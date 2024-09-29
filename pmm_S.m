function [G,W,F,H,info] = pmm_S(inputfile,opts,windowSize,proximityThreshold,q)
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
    [F,H] = readTouchstone(inputfile,opts);
    proximityThreshold=ceil(size(H,1)/2);
    if nargin < 5
        q = 2*(countPeaksAndValleys3D(H,windowSize,proximityThreshold));
        q = max(1,q);
        q = min(q,50);
        fprintf("初始阶数为%d \n",q)  
    end
    opts.q = q;
%     if size(H,1)*q<500
%         opts.q = 2*q;
%     else
%         opts.q = q;
%     end
    [F,H] = freqinterp(F,H,opts);

    opts = select_method(size(H,1),q,opts);
    check_methods(opts);
    %opts
    %归一化F
    [F1,H1,scale] = datascale(F,H,opts);
    G = []; W = [];
    if isfield(opts,'Func')
        %fprintf('\n');
        info = [];
        for k = 1:length(opts.Func)
            [G,W,info_k] = func_call(opts.Func{k},G,W,F1,H1,opts);
            info{k} = info_k;
            if ~info_k.success
               break;
            end
%             if if opts.enforceDC == 1 %原本是if opts.enforceDC == 1，这里我认为只有原本无源才能收敛，得两分
%                 G.D = real(H1(:,:,1) + G.C * (G.A \ G.B));
%             end
        end
    end

    %% Step 4: scale back
    G = ss_scale(G,scale);
    if exist('W','var') && ~isempty(W)
        W = ss_scale(W,scale);
        ss_export(G,W,subcktname,outputfile,opts);
        ss_verify(G,W,F,H,opts);
    else
        W = G;
        W.C = W.C * 0;
        ss_export(G,W,subcktname,outputfile,opts);
    end
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
        errorWithPercent = sprintf('%.5f%%', info{c}.error * 100);
        fprintf('%-15s %-15.5f %-15s %-15.5f %-15.5f %-10s\n', ...
            info{c}.func, info{c}.time, errorWithPercent, info{c}.dc_error, info{c}.k_accuracy, info{c}.passivity);
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
        dataLine = sprintf('%-15s %-15.5f %-15.5f %-15.5f %-15.5f %-10s\n', ...
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
function [G,W,info] = func_call(funcname,G,W,F1,H1,opts)
    global pmm_methods;

    func        = pmm_methods.(funcname).func;
    passive_in  = pmm_methods.(funcname).passive_in;
    passive_out = pmm_methods.(funcname).passive_out;

    if  opts.Sample == 1
        t0 = cputime;
        delta=5;
        nPoints = size(H1, 3);  
        allChangePoints = false(nPoints, 1);  

        % 计算实部零交叉点并标记变更点
        H2_real = real(diff(H1));  
        for i = 1:size(H2_real, 1)
            for j = 1:size(H2_real, 2)
                zeroCrossings = find(diff(sign(H2_real(i, j, :))));  % 实部的零交叉点
                changePoints = zeroCrossings + 1; 
                for k = 1:length(changePoints)
                    lowerBound = max(changePoints(k) - delta, 1);  % 确保下限不小于1
                    upperBound = min(changePoints(k) + delta, nPoints);  % 确保上限不大于nPoints
                    allChangePoints(lowerBound:upperBound) = true;  % 标记变化点
                end
            end
        end

        % 计算虚部零交叉点并标记变更点
        H2_imag = imag(diff(H1));  
        for i = 1:size(H2_imag, 1)
            for j = 1:size(H2_imag, 2)
                zeroCrossings = find(diff(sign(H2_imag(i, j, :))));  % 虚部的零交叉点
                changePoints = zeroCrossings + 1; 
                for k = 1:length(changePoints)
                    lowerBound = max(changePoints(k) - delta, 1);  % 确保下限不小于1
                    upperBound = min(changePoints(k) + delta, nPoints);  % 确保上限不大于nPoints
                    allChangePoints(lowerBound:upperBound) = true;  % 标记变化点
                end
            end
        end

        odd = (1:20:nPoints)';
        allChangePoints(odd) = true;  

        % 提取所有唯一的变更点索引
        uniqueChangePoints = find(allChangePoints);  
 
        H2=H1(:,:,uniqueChangePoints);
        F2=F1(uniqueChangePoints,1);


        [G,W] = func(G,W,F2,H2,opts);
        t = cputime - t0;
    else
        t0 = cputime;
        [G,W] = func(G,W,F1,H1,opts);
        t = cputime - t0;
    end
    

    [r2,f2] = passivity_violation(G);
    if passive_out && ~isempty(r2)
        fprintf('Warning: %s claims but fails to ensure passivity of the system.\n', funcname);
        fprintf('Check the following frequencies for diagnosis.\n');
        fprintf('%e\n', r2/2/pi);
        info.success = 0;
    else
        info.success = 1;
    end

    info.func = func2str(func);
    info.time = t;
    info.error = norm_error(G,F1,H1,opts);
    info.k_accuracy = info.error*opts.q; 
    if isempty(r2)
        info.passivity = 'passive';
    else
        info.passivity = 'non-passive';
    end

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
    [filepath,filename,ext] = fileparts(inputfile);
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

function [err]=norm_error(S,F,H1,opts)
    H2=ss_xf(S,F);
    N = size(H2, 3);  % 获取频率点数量
    err_numerator = 0;
    err_denominator = 0;
    for p = 1:N
        % 计算误差的分子部分: norm(F_p - S_p)
        err_numerator = err_numerator + norm(H1(:,:,p) - H2(:,:,p), 2);
        % 计算误差的分母部分: norm(S_p)
        err_denominator = err_denominator + norm(H1(:,:,p), 2);
    end
    err = err_numerator / err_denominator;
end


