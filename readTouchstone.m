function [freq,Scatter_params,param_type,impedance] = readTouchstone(filename,opts)
% readtouchstone 读取 Touchstone 文件并提取网络参数。
%   [freq, scatter_params, param_type, impedance] = readtouchstone(filename, opts)
%   读取指定的 Touchstone 格式文件，并返回频率、散射参数、参数类型和参考阻抗。
%
%   输入参数：
%       filename - Touchstone 文件的路径
%       opts     - 可选参数结构体，包含其他选项
%
%   输出参数：
%       freq         - 频率点
%       scatter_params - 提取的散射参数
%       param_type   - 网络参数类型（s、y、z、h、g）
%       impedance    - 参考阻抗值（默认 50 Ω）
    tic 
    if nargin < 2
       opts=pmm_default;
       opts.parametertype = 'S';
       opts.enforceDC = 1;
    end

    fileID = fopen(filename, 'r');

    fseek(fileID, 0, 'eof');  % 移动到文件末尾
    file_size = ftell(fileID);  % 获取文件大小
    fseek(fileID, 0, 'bof');  % 移动回文件开头

    % 初始化变量
    fre_index = 0;
    bytes_read = 0;
    param_type = 'S';  % 默认参数类型为 S 参数
    freq_scale = 1;  % 频率缩放系数，默认 Hz
    impedance = 50;  % 默认参考阻抗为 50 欧姆
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

    % 预定义 current_data 数组
    num_params = port_num^2 * 2;  % 总数据量 (port_num^2 个复数值，每个复数值有两个部分)
    current_data = zeros(1, num_params);  % 预定义数组，大小为 1 x port_num^2 * 2
    
    while ~feof(fileID)
        line = fgetl(fileID);
        
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
            continue;
        end
        
        % 提取频率和网络参数数据
        data = sscanf(line, '%f');
        if ~isempty(data)
            fre_index = fre_index+1;
            if fre_index == 1
                temp_freq = data(1) * freq_scale;
                bytes_read = bytes_read + length(line) + 1;  % 加上换行符
            end

            if fre_index ==2
                pred_num = ceil(file_size/bytes_read*1.1);
                freq = zeros(pred_num,1);
                Scatter_params = zeros(port_num,port_num,pred_num);
                freq(1,:)=temp_freq;
                Scatter_params(:,:,1)=Scatter_param;
            end
            if fre_index ~= 1
                freq(fre_index,:) = data(1) * freq_scale;  % 频率
            end
            complex_data = [];
            
             % 读取数据
            current_data(:) = 0;  % 清空数组
            current_data(1:length(data)-1) = data(2:end);  % 存储数据

            % 记录有效数据的长度
            actual_length = length(data) - 1;  % 初始有效长度
            
            % 如果当前行的数据不够，则继续读取下一行
            while actual_length < num_params
                line = fgetl(fileID);
                next_data = sscanf(line, '%f');
                if fre_index==1
                    bytes_read = bytes_read + length(line) + 1;  % 加上换行符
                end
                if ~isempty(next_data)
                    % 将新数据添加到 current_data 中，基于有效长度的位置
                    current_data(actual_length + 1 : actual_length + length(next_data)) = next_data;
                    % 更新有效数据的长度
                    actual_length = actual_length + length(next_data);  % 增加有效长度
                end
            end
 
            % 提取网络参数（根据文件中的数据格式处理复数值）
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
            if fre_index == 1
                Scatter_param = convertToS(params,param_type,impedance);
            else
                Scatter_params(:,:,fre_index) = convertToS(params,param_type,impedance);
            end
        end
    end
    fclose(fileID);

    freq = freq(1:fre_index,:);
    Scatter_params = Scatter_params(:,:,1:fre_index);
    
    % 转换为对称阵
    Scatter_params = (Scatter_params + permute(Scatter_params, [1, 2, 3])) / 2;
    %Scatter_params=Scatter_params_sym;
    % 去除高频点
    fmax = optget(opts,'fmax',1e100);
    ix=find(freq>fmax);
    if ~isempty(ix)
        freq=freq(1:ix);
        Scatter_params=Scatter_params(:,:,1:ix);
    end
    [freq, indices] = sort(freq);
    Scatter_params  = Scatter_params(:,:,indices);
    fprintf("读取文件耗时 %.8f s\n",toc) 
end

function [S_params] = convertToS(params, param_type, Z0)
    % 将其他参数类型 (Z, Y, H, G) 转换为 S 参数
    % params: 输入的网络参数矩阵
    % param_type: 参数类型 ('Z', 'Y', 'H', 'G')
    % Z0: 参考阻抗 (默认值为 50 欧姆)
    
    if nargin < 3
        Z0 = 50;  % 默认参考阻抗
    end
    
    % 计算参考导纳
    Y0 = 1 / Z0;
    
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
