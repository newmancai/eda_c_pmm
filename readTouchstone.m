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

    if nargin < 2
       opts=pmm_default;
       opts.parametertype = 'S';
       opts.enforceDC = 1;
    end

    fileID = fopen(filename, 'r');
    % 初始化变量
    freq = [];
    Scatter_params = [];
    port_num = 0;
    param_type = 'S';  % 默认参数类型为 S 参数
    freq_scale = 1;  % 频率缩放系数，默认 Hz
    impedance = 50;  % 默认参考阻抗为 50 欧姆
    data_format = 'RI';  % 默认数据格式为 RI (Real + Imaginary)

    % 从文件后缀 snp 中提取端口数 n
    [~, ~, ext] = fileparts(filename);
    if startsWith(ext, '.s') && length(ext) > 2
        port_num = str2double(ext(3:end-1));  % 端口数
    else
        error('文件后缀无效，无法确定端口数量');
    end
    
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

            freq = [freq; data(1) * freq_scale];  % 频率

            complex_data = [];
            
            % 需要读取的总数据量 (port_num^2 个复数值，每个复数值有两个部分)
            num_params = port_num^2 * 2;
            current_data = data(2:end);  % 去掉频率值
            
            % 如果当前行的数据不够，则继续读取下一行
            while length(current_data) < num_params
                line = fgetl(fileID);
                next_data = sscanf(line, '%f');
                current_data = [current_data; next_data];
            end
            
            % 取出足够的复数数据
            current_data = current_data(1:num_params);
  
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
            Scatter_param = convertToS(params,param_type,impedance);
            Scatter_params = [Scatter_params; Scatter_param];
        end
    end
    
    fclose(fileID);
    
    % 将网络参数重塑为 (port_num, port_num, frequency_points) 的三维矩阵
    num_freq = length(freq);
    Scatter_params = permute(reshape(Scatter_params, [port_num, num_freq, port_num]), [1, 3, 2]);
    % 转换为对称阵
    ns=size(Scatter_params,3);
    for c=1:ns
          Scatter_params_sym(:,:,c)=(Scatter_params(:,:,c)+Scatter_params(:,:,c).')/2;
    end
    Scatter_params=Scatter_params_sym;
    % 去除高频点
    fmax = optget(opts,'fmax',1e100);
    ix=find(freq>fmax);
    if ~isempty(ix)
        freq=freq(1:ix);
        Scatter_params=Scatter_params(:,:,1:ix);
    end
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
