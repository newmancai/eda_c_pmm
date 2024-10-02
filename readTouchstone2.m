function [freq, Scatter_params, param_type, impedance] = readTouchstone2(filename, opts)
    % readtouchstone 读取 Touchstone 文件并提取网络参数。
    %   [freq, scatter_params, param_type, impedance] = readtouchstone(filename, opts)
    %   读取指定的 Touchstone 格式文件，并返回频率、散射参数、参数类型和参考阻抗。

    tic;  % 开始计时
    if nargin < 2
        opts = pmm_default;
        opts.parametertype = 'S';
        opts.enforceDC = 1;
    end

    fileID = fopen(filename, 'r');
    raw_data = fread(fileID, '*char')';  
    fclose(fileID);

    % 初始化变量
    line_idx = 1;
    param_type = 'S';  % 默认参数类型为 S 参数
    freq_scale = 1;    % 频率缩放系数，默认 Hz
    impedance = 50;    % 默认参考阻抗为 50 欧姆
    data_format = 'RI';  % 默认数据格式为 RI (Real + Imaginary)

   % 从文件后缀 snp 中提取端口数 n
    [~, ~, ext] = fileparts(filename);
    if startsWith(ext, '.s') && length(ext) > 2
        port_num = str2double(ext(3:end-1));  % 端口数
    elseif startsWith(ext, '.S') && length(ext) > 2
        port_num = str2double(ext(3:end-1));  % 端口数
    else
        error('文件后缀无效，无法确定端口数量');
    end

    % 初始化num_params
    num_params = port_num^2 * 2;  % 总数据量 (port_num^2 个复数值，每个复数值有两个部分)

    % 找到文件头部和注释部分的结束行（跳过注释和头部信息）
    length_raw_data = length(raw_data);
    for i = 1:length_raw_data
        if raw_data(i) == newline  % 检查是否为换行符
            line = raw_data(line_idx:i-1);  % 截取当前行的字符串
            line_idx = i + 1;  % 增加行数

            % 跳过注释行（以!号开头的行）
            if startsWith(line, '!')
                continue;
            end

            % 处理文件头部信息，获取频率单位、参数类型和参考阻抗
            if startsWith(line, '#')
                parts = split(line);

                % 解析频率单位
                freq_unit = lower(parts{2});  % 频率单位 (Hz, MHz, GHz)
                if strcmp(freq_unit, 'ghz')
                    freq_scale = 1e9;
                elseif strcmp(freq_unit, 'mhz')
                    freq_scale = 1e6;
                elseif strcmp(freq_unit, 'khz')
                    freq_scale = 1e3;
                elseif strcmp(freq_unit, 'hz')
                    freq_scale = 1;
                end

                % 解析参数类型 (S, Y, Z, H, G)
                param_type = upper(parts{3});  % 网络参数类型（S, Y, Z, H, G）

                % 解析复数数据的表示形式 (RI, MA, DB)
                data_format = upper(parts{4});

                % 解析参考阻抗
                if length(parts) >= 6
                    impedance = str2double(parts{6});  % 参考阻抗
                end
                new_line_index = i;
                break;
            end
        end
    end

    start_index = new_line_index + 1;  % 数据行开始的位置
    % 读取文件后，提取数据行
    data_threshold = 5e8; % 设置一个阈值大小，超过该值时进行分块读取
    
    if length_raw_data > data_threshold
        % 文件特别大，分块读取
        block_size = 5e5;  % 每次处理1亿个字符的块
        full_data = [];
        current_pos = start_index;
        while current_pos < length_raw_data
            % 计算预期的当前块结束位置
            expected_end_pos = min(current_pos + block_size - 1,length_raw_data);

            % 寻找最近的换行符，确保不截断数据
            newline_pos = find(raw_data(current_pos:expected_end_pos) == newline, 1, 'last');

            if isempty(newline_pos)
                % 没有找到换行符，可能块刚好结束在换行符处
                end_pos = expected_end_pos;
            else
                % 找到换行符，调整end_pos
                end_pos = current_pos + newline_pos - 1;
            end

            % 按行拆分数据并预分配有效行的空间
            block_data = raw_data(current_pos:end_pos);
            if contains(block_data, '!')
                lines = strsplit(block_data,'\n');
                cleaned_block = repmat("", length(lines), 1);  % 预分配空间
                cleaned_index = 0;  % 用于跟踪有效的非注释行索引

                for i = 1:length(lines)
                    % 删除以!开头的注释行
                    if ~startsWith(lines{i}, '!')
                        cleaned_index = cleaned_index + 1;  % 增加有效行计数
                        cleaned_block(cleaned_index) = lines{i};  % 仅保留非注释行
                    end
                end

                % 删除未使用的行
                cleaned_block = cleaned_block(1:cleaned_index);

                % 将清理后的块合并为字符串
                cleaned_data = strjoin(cleaned_block, newline);

                % 使用 sscanf 直接处理合并后的字符串
                block_data = sscanf(cleaned_data, '%f');
            else
                block_data = sscanf(block_data, '%f');
            end
            %fprintf("Size of block_data: %d \n", size(block_data, 1));

            % 累加到 full_data 中
            full_data = [full_data; block_data];
            % 更新当前位置
            current_pos = end_pos + 1;
        end
        clear raw_data;
    else
        % 文件不大，直接处理
        raw_data = raw_data(start_index:end);  % 提取数据行
        % 按行拆分数据并预分配有效行的空间
        lines = strsplit(raw_data, newline);
        clear raw_data;
        cleaned_block = repmat("", length(lines), 1);  % 预分配空间
        cleaned_index = 0;  % 用于跟踪有效的非注释行索引

        for i = 1:length(lines)
            % 删除以!开头的注释行
            if ~startsWith(lines{i}, '!')
                cleaned_index = cleaned_index + 1;  % 增加有效行计数
                cleaned_block(cleaned_index) = lines{i};  % 仅保留非注释行
            end
        end

        % 删除未使用的行
        cleaned_block = cleaned_block(1:cleaned_index);

        % 将清理后的块合并为字符串
        cleaned_data = strjoin(cleaned_block, newline);
        
        % 使用 sscanf 直接处理合并后的字符串
        full_data = sscanf(cleaned_data, '%f');
    end

    % 现在你已经有了所有的数据，可以批量处理 full_data
    data_count = length(full_data);
    num_points = floor(data_count / (num_params + 1));  % 计算总频率点的数量

     % 初始化输出变量
    freq = zeros(num_points, 1);
    Scatter_params = zeros(port_num, port_num, num_points);
    
    % 解析频率和参数数据
    data_idx = 1;
    for i = 1:num_points
        freq(i) = full_data(data_idx) * freq_scale;
        data_idx = data_idx + 1;
        current_data = full_data(data_idx:data_idx + num_params - 1);
        data_idx = data_idx + num_params;

        % 处理复数数据
        switch data_format
            case 'RI'  % 实部 + 虚部
                real_part = current_data(1:2:end-1);
                imag_part = current_data(2:2:end);
                complex_data = real_part + 1i * imag_part;
            case 'MA'  % 幅度 + 相位（角度）
                magnitude = current_data(1:2:end-1);
                angle = deg2rad(current_data(2:2:end));  
                complex_data = magnitude .* exp(1i * angle);
            case 'DB'  % dB + 相位（角度）
                dB_value = current_data(1:2:end-1);
                angle = deg2rad(current_data(2:2:end));  
                magnitude = 10.^(dB_value / 20);  
                complex_data = magnitude .* exp(1i * angle);
        end

        % 将提取的数据重塑为 (port_num x port_num) 矩阵
        params = reshape(complex_data, [port_num, port_num]);
        Scatter_params(:,:,i) = convertToS(params, param_type, impedance);
    end
       
    % 转换为对称阵
    Scatter_params = (Scatter_params + permute(Scatter_params, [1, 2, 3])) / 2;
    
    % 去除高频点
    fmax = optget(opts,'fmax',1e100);
    ix = find(freq > fmax);
    if ~isempty(ix)
        freq = freq(1:ix);
        Scatter_params = Scatter_params(:,:,1:ix);
    end

    [freq, indices] = sort(freq);
    Scatter_params  = Scatter_params(:,:,indices);
    
    fprintf("读取文件耗时 %.8f s\n", toc);
end

function [S_params] = convertToS(params, param_type, Z0)
    % 将其他参数类型 (Z, Y, H, G) 转换为 S 参数
    if nargin < 3
        Z0 = 50;  % 默认参考阻抗
    end
    
    Y0 = 1 / Z0;  % 参考导纳
    
    % 参数类型转换
    switch param_type
        case 'S'
            S_params = params;
        case 'Z'
            S_params = z2s(params,Z0);
        case 'Y'
            S_params = y2s(params,Z0);
        case 'H'
            S_params = h2s(params,Z0);
        case 'G'
            S_params = g2h(params);
            S_params = h2s(S_params,Z0);
        otherwise
            error('不支持的参数类型');
    end
end
