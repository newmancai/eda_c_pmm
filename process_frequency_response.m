function real_num = process_frequency_response(y11_si, s, poles, residues, d,all_type)
    % PROCESS_FREQUENCY_RESPONSE - 去除常数项、一阶项和实极点的贡献，添加负频率对称性并进行 SVD 分析
    %
    % 输入参数：
    % y11_si    - 原始频率响应样本，大小为 MxMxs
    % s         - 频率样本点向量
    % poles     - 极点矩阵，大小为 (Nq, 1)，包含 Nq 个复数极点
    % residues  - 留数矩阵，大小为 MxMxNq，M 为维度，Nq 为极点数量
    % d         - 常数项，大小为 MxM
    %
    % 输出参数：
    % real_num  - 实极点的数量和 SVD 阶数分析后的修正值
    
    tic
    [M, ~, Nf] = size(y11_si);  
    [Nq, ~] = size(poles);      
    
    y11_si_hat = y11_si; 
    real_num = 0;
    
    % 第一步：去除常数项 d
    for i = 1:Nf
        y11_si_hat(:,:,i) = y11_si_hat(:,:,i) - d; 
    end

    % 第二步：去除实极点的贡献
    for k = 1:Nq
        if isreal(poles(k))  
            real_num = real_num + 1;
            for i = 1:Nf
                y11_si_hat(:,:,i) = y11_si_hat(:,:,i) - residues(:,:,k) ./ (s(i) - poles(k));
            end
        end
    end

    % 第五步：添加负频率的共轭对称性
    N = Nf;  
    y11_neg_freq = conj(y11_si_hat(:,:,end-1:-1:1));  
    y11_full_spectrum = cat(3, y11_si_hat, y11_neg_freq);  

    % 对每个 MxM 的频率响应矩阵进行 IFFT
    if all_type ==1
        x = zeros(M, M, 2*N-1);
        for i = 1:M
            for j = 1:M
                spectrum_ij = squeeze(y11_full_spectrum(i,j,:));
                x(i,j,:) = ifft(spectrum_ij, 'symmetric');
            end
        end
    else
        x = zeros(1, 1, 2*N-1);
        for i = 1:1
            for j = 1:1
                spectrum_ij = squeeze(y11_full_spectrum(i,j,:));
                x(i,j,:) = ifft(spectrum_ij, 'symmetric');
            end
        end
    end
    
    % 第六步：构建系数矩阵 R 并进行 SVD 分析
    x_1d = squeeze(x(1, 1, :));  
    p = 40;  
    R = build_matrix(x_1d, p);  
    S = svds(R, p);
    min_threshold = 0.00000001;  % 阈值范围的最小值
    max_threshold = 0.01;        % 阈值范围的最大值
    for i = 2:2:p
        if S(i)/S(1) < calculate_threshold(real_num, 1, min_threshold, max_threshold)  
            break;
        end
    end
    real_num = real_num + i;
    toc
end

function R = build_matrix(x_1d, p)
    % BUILD_MATRIX - 构建矩阵，用于后续 SVD 分析
    %
    % 输入参数：
    % x_1d - 一维数组，长度为 N
    % p    - 矩阵的行数
    %
    % 输出参数：
    % R - 构建的矩阵，大小为 p x (N-p-1)
    
    N = length(x_1d);  
    
    R = zeros(p, N-p-1);
    
    for i = 1:p
        R(i, :) = x_1d(i:N-p+i-2);
    end
end

function threshold = calculate_threshold(real_num, sum_num, min_threshold, max_threshold)
    % 对阈值范围取对数
    log_threshold_min = log10(min_threshold);
    log_threshold_max = log10(max_threshold);
    
    % 计算对应的阈值
    log_threshold = log_threshold_min + (log_threshold_max - log_threshold_min) * (real_num / sum_num);
    threshold = 10^log_threshold;
end