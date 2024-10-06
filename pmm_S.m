function [G,W,F,H,info] = pmm_S(inputfile,opts,windowSize,proximityThreshold)
% PMM Passive Macro Modeling (PMM) function.
%
%   [G, W, F, H, INFO] = PMM(INPUTFILE, Q, OPTS) performs passive macro modeling
%   on input data. It begins by loading data from the input file, performing 
%   frequency interpolation, selecting appropriate modeling methods, and scaling 
%   the data before final processing.
%
%   Inputs:
%       INPUTFILE - The input file containing the system data.
%       Q - The number of poles used in the model.
%       OPTS - A structure containing options for the PMM process.
%
%   Outputs:
%       G - State-space model.
%       W - Frequency response.
%       F - Frequency vector.
%       H - Transfer function matrix.
%       INFO - Additional details about the processing steps.
%
%   Example:
%       [G, W, F, H, info] = pmm('symind.s2p', 5, opts);
%
%   See also: LDSTONE, FREQINTERP, DATASCALE
%
%   Author: Cai Yuyang
%   Email: [email-redacted]
%   Modified: 23-Sep-2024 - Initial version.
% Copyright 2012, Zuochang Ye, zuochang@tsinghua.edu.cn
% --- Function code starts here ---
    if nargin < 3
        opts = pmm_default;
    end
    %opts.q = q;

    [inputfile, outputfile, subcktname] = ParseInput(inputfile);

    %% Step 0: load data
    [F,H] = readTouchstone2(inputfile,opts);
    if nargin < 4
        proximityThreshold = 30;%合并阈值设为和端口数量相关，不然大电路无法运行
    end
    
%     if nargin < 5
    [valleypeakIndices,q] = countPeaksAndValleys3D(H,windowSize,proximityThreshold);
%     if size(H,1)>60
%         opts.jk = 1;
%     end
    q = min(max(ceil(2*opts.jk*q), 6), 50);
    q = q + mod(q, 2);
    opts.q = q;
    fprintf("初始阶数为 %d \n",q);
    [F,H] = freqinterp(F,H,opts);
    port_num = size(H,1);
    opts = select_method(size(H,1),q,opts);
    check_methods(opts);
%     opts.q = opts.q*2;
    
    %归一化F
    [F1,H1,scale] = datascale(F,H,opts);
    

    G = []; W = [];
    min_k_accuracy = [];
    allG = [];
    count_vf = 0;
    if isfield(opts,'Func')
        fprintf('\n');
        info = [];
        newFunc = opts.Func;  
        k = 1; 
        threshold = 1e-5;
        while k <= length(newFunc)
            [G_k,~,info_k,F2,H2] = func_call(newFunc{k},G,W,F1,H1,opts,valleypeakIndices,port_num,threshold);
            allG{k}= G_k;
            info{k} = info_k;
            if ~info_k.success
               break;
            end
            G = G_k;
            [residues,poles]=ss2pr(G.A,G.B,G.C);
            Hinf=G.D;
            fprintf("目前的error: %d \n",info{k}.error);
            % 计算推荐阶数并运行        
            min_k_accuracy = [min_k_accuracy,info{k}.k_accuracy];
            order = process_frequency_response(H1,F1*2*pi*1j, poles, residues,Hinf);
            fprintf("推荐阶数：%d \n",order);
            if (info{k}.error > 0.1)||(order/size(poles,1)>1.1)
                count_vf = count_vf + 1;
                if (count_vf > 10)%看看结果，如何优化？
                    break;
                end
                if(info{k}.error > 0.1)
                    if(info{k}.error > 0.2)
                        order = max(size(poles,1)+ 20,order) ;
                    else
                        order = max(size(poles,1)+ 10,order) ;
                    end    
                    threshold = threshold / 10;
                end
                if order>(size(F2,1)-2)
                    opts.sample_add = 1;
                end
                newFunc = [newFunc(1), {'VF'}, newFunc(2:end)];
                opts.q = order;
            else
                if size(H1,1)^2*opts.q <65536
                    G = optimizeSystem(G, F2, H2, opts);
                    info_k.error = norm_error(G,F1,H1);
                    info_k.k_accuracy = info_k.error*opts.q;
                    H0 = H1(:,:,1);
                    HH = G.D - G.C * (G.A \ G.B);
                    info_k.dc_error = max(vec(abs(H0 - HH)));
                end
            end
            opts.Func = newFunc;
            %if opts.enforceDC == 1
            %     G.D = H1(:,:,1) + G.C * (G.A \ G.B);
            %end
            k = k + 1;
        end
    end

    [~, minIndex] = min(min_k_accuracy);
    G = allG{minIndex};

    %% Step 4: scale back
    G = ss_scale(G,scale);
    if exist('W','var') && ~isempty(W)
        W = ss_scale(W,scale);
        %ss_export(G,W,subcktname,outputfile,opts);
        ss_verify(G,W,F,H,opts);
    else
        W = G;
        W.C = W.C * 0;
        %ss_export(G,W,subcktname,outputfile,opts);
    end
    print_info_additional(info,extract_filename(inputfile));
    %print_info(info)
    %print_info_to_file(info,inputfile)
end

%{
    @brief Function to print detailed information about the PMM process.
    
    This function prints the function names, processing time, errors, and passivity status
    of each step during the PMM execution.
    
    @param info A structure containing detailed information about each function used 
    in the PMM process.
%}
function print_info(info)
%sss
    col_widths = [15, 15, 15, 15, 15, 10]; % 各列宽度
    total_length = sum(col_widths) + length(col_widths) - 1; % 计算总长度（包含空格）
    fprintf('%-15s %-15s %-15s %-15s %-15s %-10s\n', ...
        'FuncName', 'Time', 'Error(Norm)', 'Error(DC)', 'K_accuracy', 'Passivity');
    fprintf('%s\n', repmat('-', 1, total_length));
    
    % 输出每个信息的内容
    for c = 1:length(info)
        fprintf('%-15s %-15d %-15d %-15d %-15d %-10s\n', ...
            info{c}.func, info{c}.time, info{c}.error, info{c}.dc_error, info{c}.k_accuracy, info{c}.passivity);
    end
end

function print_info_additional(info,filename)
%sss
    col_widths = [15 ,10, 10, 15, 20, 15, 20, 10]; % 各列宽度
    total_length = sum(col_widths) + length(col_widths) - 1; % 计算总长度（包含空格）
    fprintf('%-15s %-10s %-10s %-15s %-20s %-20s %-15s %-10s\n', ...
        'DocName','FuncName','Order','Time', 'Error(Norm)', 'Error(DC)', 'K_accuracy', 'Passivity');
    fprintf('%s\n', repmat('-', 1, total_length));
    
    % 输出每个信息的内容
    for c = 1:length(info)
        errorWithPercent = sprintf('%.9f%%', info{c}.error * 100);
        fprintf('%-15s %-10s %-10d %-15.5f %-20s %-20.11f %-15.8f %-10s\n', ...
            filename,info{c}.func,info{c}.order,info{c}.time, errorWithPercent, info{c}.dc_error, info{c}.k_accuracy, info{c}.passivity);
    end
end

function print_info_to_file(info, input_filename)
% PRINT_INFO_TO_FILE Outputs structured information to a log file.
%   PRINT_INFO_TO_FILE(INFO, INPUT_FILENAME) takes a cell array of structures
%   INFO containing fields such as func, time, error, dc_error, k_accuracy,
%   and passivity, and an INPUT_FILENAME to derive the output filename.
%   The output filename will include a timestamp and have a .log extension.
%
%   The resulting log file will contain a header and the corresponding
%   information formatted in aligned columns.

    [~, name, ~] = fileparts(input_filename);
    timestamp = datestr(now, 'dd_HHMMSS');
    output_filename = [name, '_', timestamp, '.log'];
    fid = fopen(output_filename, 'w');

    if fid == -1
        error('无法打开文件 %s', output_filename);
    end
    
    col_widths = [15, 15, 15, 15, 15, 10];
    total_length = sum(col_widths) + length(col_widths) - 1;

    header = sprintf('%-15s %-15s %-15s %-15s %-15s %-10s\n', ...
        'FuncName', 'Time', 'Error(Norm)', 'Error(DC)', 'K_accuracy', 'Passivity');
    fprintf(fid, header);
    fprintf(fid, '%s\n', repmat('-', 1, total_length));
    
    for c = 1:length(info)
        dataLine = sprintf('%-15s %-15e %-15e %-15e %-15e %-10s\n', ...
            info{c}.func, info{c}.time, info{c}.error, info{c}.dc_error, info{c}.k_accuracy, info{c}.passivity);
        fprintf(fid, dataLine);
    end
    
    fclose(fid);
end


%{
    @brief Calls specific PMM methods to process the data.
    
    This function calls the appropriate method (specified by the user) to process 
    the input data and update the state-space matrices G and W.
    
    @param funcname The name of the method to be called.
    @param G State-space model matrix G.
    @param W Frequency response matrix W.
    @param F1 Interpolated frequency data.
    @param H1 Interpolated transfer function data.
    @param opts Options structure containing parameters for the method.
    
    @return [G, W, info] Updated state-space model and frequency response.
%}


function [H2, F2, uniqueChangePoints, valleypeak_sorted] = processH1(H1, F1, valleypeakIndices, opts)
    % 输入：
    % H1 - 输入的三维矩阵
    % F1 - 输入的频率数组
    % valleypeakIndices - 输入的峰值索引
    % opts.Sample - 采样选项（如果为1，进行采样）

    % Step 1: 计算 valleypeak 和 H1_choose
    valleypeak = 2 .* pi .* F1(valleypeakIndices);
    H1_choose = sum(sum(H1(:,:,valleypeakIndices)));
    
    % Step 2: 按绝对值降序排列
    [H1_choose_sorted, indices] = sort(abs(H1_choose(:)), 'descend');
    valleypeak_sorted = valleypeak(indices);
    
    % Step 3: 判断 opts.Sample 的值，如果为1进行采样
    if opts.Sample == 1
        delta = 10;
        if opts.sample_add 
            delta = delta +5;
        end
        delta2 = 2;
        delta3 = 50;
        allChangePoints = [];

        % Step 4: 计算所有的变化点
        for i = 1:size(valleypeakIndices, 1)
            lowerBound = max(valleypeakIndices(i) - delta, 1); % 确保不小于0
            upperBound = min(valleypeakIndices(i) + delta, size(H1, 3)); % 确保不大于H1的第三维度
            allChangePoints = [allChangePoints; (lowerBound:delta2:upperBound)'];
        end

        % Step 5: 添加额外的变化点
        odd = (1:delta3:size(H1, 3))';
        allChangePoints = [allChangePoints; odd];

        % Step 6: 找到所有唯一的变化点
        uniqueChangePoints = unique(allChangePoints);

        % Step 7: 使用变化点提取 H2 和 F2
        H2 = H1(:,:,uniqueChangePoints);
        F2 = F1(uniqueChangePoints, 1);
    else
        % 如果 opts.Sample 不为1，返回空矩阵
        H2 = [];
        F2 = [];
        uniqueChangePoints = [];
    end
end

function [G,W,info,F2,H2] = func_call(funcname,G,W,F1,H1,opts,valleypeakIndices,port_num,threshold)
    global pmm_methods;

    func        = pmm_methods.(funcname).func;
    passive_in  = pmm_methods.(funcname).passive_in;
    passive_out = pmm_methods.(funcname).passive_out;
    
    
    valleypeak=2.*pi.*F1(valleypeakIndices);
    H1_choose = sum(sum(H1(:,:,valleypeakIndices)));
    [H1_choose_sorted, indices] = sort(abs(H1_choose(:)), 'descend');
    valleypeak = valleypeak(indices);
    
    if  opts.Sample == 1
        delta=10;
        if opts.sample_add 
            delta = delta +5;
        end
        delta2=2;
        delta3=50;
        allChangePoints = [] ;

        for i =1:size(valleypeakIndices,1)
            lowerBound = max(valleypeakIndices(i) - delta, 1); % 确保不小于0
            upperBound = min(valleypeakIndices(i) + delta, size(H1,3)); % 确保不大于endValue
            allChangePoints = [allChangePoints; (lowerBound:delta2:upperBound)'];
        end

        odd=(1:delta3:size(H1,3))';
        allChangePoints = [allChangePoints;odd];
        uniqueChangePoints = unique(allChangePoints);

        H2=H1(:,:,uniqueChangePoints);
        F2=F1(uniqueChangePoints,1);
        opts.q = min(size(F2,1)-2,opts.q); %防止QR分解后的下三角为空
        if size(H2,3)>=6
            Nc = size(H2, 1);
            reshaped_H2 = reshape(H2, Nc * Nc, size(H2, 3))';
            concatenated_H2 = [real(reshaped_H2); imag(reshaped_H2)];

            rho = 2 * Nc;
            [U_truncated, S_truncated, V_truncated] = svds(concatenated_H2, rho);

            singular_values = diag(S_truncated);
            effective_indices = find(singular_values > threshold);
            fprintf("列空间降维为：%d 维 \n",size(effective_indices,1));

            U_effective = U_truncated(:, effective_indices);
            S_effective = S_truncated(effective_indices, effective_indices);
            V_effective = V_truncated(:, effective_indices);

            identity_matrix = eye(size(H2, 3));
            complex_identity = [identity_matrix, 1i * identity_matrix];

            W_temp_effective = U_effective * S_effective;
            W_output = complex_identity * W_temp_effective;
            opts.enableSVD = 1;
        end
        %X = W * V_rho'; %低秩下成立，现在拟合列空间的基向量


        %     size(H2,3)
        %     size(H1,3)
        t0 = clock;
        if(strcmp(func2str(func),'VF'))
            if opts.enableSVD
                [G,W] = func(G,W,F2,W_output,opts,valleypeak);
                I_P = eye(port_num);
                % Step 2: 计算 ABCD 矩阵
                G.A = kron(I_P,G.A);  % A = I_P ⊗ A_w
                G.B = kron(I_P,G.B);  % B = I_P ⊗ b_w
                C = zeros(port_num,size(G.C,2));
                D = zeros(port_num,port_num);
                n_c = size(G.C,2);
                for j = 1:port_num
                    V_j = V_effective((j-1)*port_num + 1:j*port_num, :);  % 取出 V 的 P 行
                    C(:,(j-1)*n_c + 1:j*n_c) = V_j*G.C;
                    D(:,j)=V_j * G.D;
                end
                G.C = C;
                G.D = D;
            else
                [G,W] = func(G,W,F2,H2,opts,valleypeak);
            end
        else
            [G,W] = func(G,W,F2,H2,opts);
        end
        t = etime(clock,t0);
    end
    
    if  opts.Sample == 2
        t0 = clock;
        
        [G,W] = func(G,W,F1,H1,opts,valleypeak);
        t = etime(clock,t0);
    end
    
     
%     [r2,~] = passivity_violation(G);
%     if passive_out && ~isempty(r2)
%         fprintf('Warning: %s claims but fails to ensure passivity of the system.\n', funcname);
%         fprintf('Check the following frequencies for diagnosis.\n');
%         fprintf('%e\n', r2/2/pi);
%         info.success = 0;
%     else
%         info.success = 1;
%     end
    info.success = 1;
    r2 = 1;
    info.func = func2str(func);
    info.time = t;
    info.error = norm_error(G,F1,H1);
    info.order = opts.q;
    info.k_accuracy = info.error*opts.q; 
    if isempty(r2)
        info.passivity = 'passive';
    else
        info.passivity = 'non-passive';
    end
   
%     G.D=H1(:,:,1)+G.C * (G.A \ G.B);
    
    H0 = H1(:,:,1);
    HH = G.D - G.C * (G.A \ G.B);
    
    info.dc_error = max(vec(abs(H0 - HH)));
end

%{
    @brief Parses the input file and sets the output file and subcircuit name.
    
    This function extracts the path, filename, and extension of the input file, and
    generates the output file and subcircuit names based on the input filename.
    
    @param inputfile The path to the input file.
    
    @return [inputfile, outputfile, subcktname] Processed input, output file path, 
    and subcircuit name.
%}
function [inputfile,outputfile,subcktname] = ParseInput(inputfile)
    [filepath,filename,~] = fileparts(inputfile);
    if isempty(filepath)
        filepath = '.';
    end
    outputfile = sprintf('%s/%s_ncss.scs', filepath, filename);
    subcktname = sprintf('%s_ncss', filename);
end
%{
    @brief Selects the appropriate method for the PMM process based on input size and options.
    
    This function selects the correct modeling method (e.g., VF, SDP, LC) based on the 
    number of variables in the system and the user-defined options.
    
    @param m Number of rows in the transfer function matrix H.
    @param q Number of poles used in the model.
    @param opts Options structure containing user settings for the method selection.
    
    @return opts Updated options structure with the selected method.
%}
function opts = select_method(m, q, opts)
    method = optget(opts, 'method', 'auto');
    if strcmp(method, 'user_defined')
        return;
    end
    opts.Func = {};
    opts.Func{1} = 'VF';
    
    if strcmp(method, 'auto')
        n = m * q;
        if n * m < 100
            opts.Func{2} = 'SDP';
        else        
            opts.Func{2} = 'LC';
            opts.Func{3} = 'DAO';        
        end
        return;
    end

    if strcmp(method, 'sdp')
        opts.Func{2} = 'SDP';
        return;
    end

    if strcmp(method, 'vf_only')
        return;
    end

    if strcmp(method, 'lc_only')
        opts.Func{2} = 'LC';    
        return;
    end

    if strcmp(method, 'asym_only')
        opts.Func{2} = 'ASYM';    
        return;
    end

    if strcmp(method, 'epm_only')
        opts.Func{2} = 'EPM';    
        return;
    end

    if strcmp(method, 'frp_only')
        opts.Func{2} = 'FRP';    
        return;
    end

    if strcmp(method, 'lc_dao')
        opts.Func{2} = 'LC';
        opts.Func{3} = 'DAO';
        return;
    end

    if strcmp(method, 'epm_dao')
        opts.Func{2} = 'EPM';    
        opts.Func{3} = 'DAO';
        return;
    end

    if strcmp(method, 'frp_dao')
        opts.Func{2} = 'FRP'; 
        opts.Func{3} = 'DAO';
        return;
    end

    if strcmp(method, 'epm2_only')
        opts.Func{2} = 'EPM2';    
        return;
    end

    if strcmp(method, 'epm3_only')
        opts.Func{2} = 'EPM3';    
        return;
    end

    if strcmp(method, 'epm4_only')
        opts.Func{2} = 'EPM4';    
        return;
    end
end

%{
    @brief Checks the validity of the methods selected in the options.
    
    This function ensures that all methods specified in the options are available 
    and installed in the system.
    
    @param opts Options structure containing user-defined methods.
    
    @throws Error if a specified method is not installed or available.
%}
function check_methods(opts)
    global pmm_methods;

    for c = 1:length(opts.Func)
        if ~isfield(pmm_methods, opts.Func{c})
            error(sprintf('Method %s not installed.', opts.Func{c}));
        end
    end
end


function new_filename = extract_filename(original_filename)
    % 提取文件名前部分并保留后缀
    % 输入:
    %   original_filename - 原始文件名 (例如 'sp125_uniform.s64p')
    % 输出:
    %   new_filename - 提取后的文件名 (例如 'sp125.s64p')

    % 使用 fileparts 分离文件名和文件后缀
    [~, name, ext] = fileparts(original_filename);

    % 找到下划线的位置
    underscore_idx = strfind(name, '_');

    % 如果找到下划线，提取下划线前面的部分，并加上文件后缀
    if ~isempty(underscore_idx)
        new_filename = [name(1:underscore_idx(1)-1), ext];
    else
        new_filename = [name, ext]; % 没有下划线的情况，直接使用原名和后缀
    end
end

function G = optimizeSystem(G, F, H, opts)
    % Input:
    % G - System structure containing matrices A, B, C, D
    % F - Input matrix for preprocessing
    % H - Output matrix for preprocessing
    % opts - Option structure

    % Check if F(1) is 0 and enforceDC option is enabled
    if F(1) == 0 && optget(opts, 'enforceDC', 0)
        [n, m] = size(G.B);  % Get the dimensions of G.B

        % Generate sparse matrix using Kronecker product
        C_sparse = kron(-G.B' * inv(G.A)', eye(m));

        % Construct equality constraints for the optimization problem
        Aeq = [C_sparse, eye(m^2)];
        H0 = H(:, :, 1);
        beq = vec(H0);

        % Remove the first elements from F and H matrices
        F = F(2:end);
        H = H(:, :, 2:end);

        % Preprocess the inputs
        [R, QG, ~] = preprocess(F, H, G, opts);

        % Set up the quadratic programming problem
        Q = R' * R;
        f = -QG' * R;
        vecG = [vec(G.C); vec(G.D)];
        % samll num
        small_value = 1e-8;  % 可以根据需要调整这个值

        % replace 0
        vecG(vecG == 0) = small_value;
        options = optimset('LargeScale', 'off', 'TypicalX', vecG, 'Display', 'off');

        % Solve the quadratic programming problem
        y = quadprog(Q, f, [], [], Aeq, beq, [], [], [], options);

        % Take the real part of the result
        y = real(y);

        % Reshape the solution into G.C and G.D matrices
        G.C = reshape(y(1:m * n), m, n);
        G.D = reshape(y((m * n + 1):(n * m + m^2)), m, m);

        % Update G.D with the current values of G.C, G.A, and G.B
        G.D = H0 + G.C * (G.A \ G.B);
    end
end
