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
    if size(H,1)>60
        opts.jk = 1;
    end
    q = min(max(ceil(2*opts.jk*q), 6), 50);
    q = q + mod(q, 2);
    opts.q = q;
    fprintf("初始阶数为 %d \n",q);
    [F,H] = freqinterp(F,H,opts);
    
    opts = select_method(size(H,1),q,opts);
    check_methods(opts);
%     opts.q = opts.q*2;
    
    %归一化F
    [F1,H1,scale] = datascale(F,H,opts);
    G = []; W = [];
    min_k_accuracy = [];
    allG = [];
    if isfield(opts,'Func')
        fprintf('\n');
        info = [];
        newFunc = opts.Func;  
        k = 1;  
        while k <= length(newFunc)
            [G_k,~,info_k] = func_call(newFunc{k},G,W,F1,H1,opts,valleypeakIndices);
            allG{k}= G_k;
            info{k} = info_k;
            if ~info_k.success
               break;
            end
            G = G_k;
            [residues,poles]=ss2pr(G.A,G.B,G.C);
            Hinf=G.D;
            % 计算推荐阶数并运行        
            min_k_accuracy = [min_k_accuracy,info{k}.k_accuracy];
            order = process_frequency_response(H1,F1*2*pi*1j, poles, residues,Hinf);
            fprintf("推荐阶数：%d \n",order);
            if (info{k}.k_accuracy/info{k}.order > 0.1)||((order*size(H1,1)^2<7e4)&&(order>size(poles,1)+2))
                if(order<size(poles,1))
                    order = size(poles,1) + 4;
                end
                newFunc = [newFunc(1), {'VF'}, newFunc(2:end)];
                opts.q = order;
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
        ss_export(G,W,subcktname,outputfile,opts);
        ss_verify(G,W,F,H,opts);
    else
        W = G;
        W.C = W.C * 0;
        ss_export(G,W,subcktname,outputfile,opts);
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
function [G,W,info] = func_call(funcname,G,W,F1,H1,opts,valleypeakIndices)
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
    delta2=2;
    delta3=50;
%     H2=real(diff(H1));
    allChangePoints = [] ;
%     for i =1:size(H2,1)
%         for j =1:size(H2,2)
%             zeroCrossings = find(diff(sign(H2(i,j,:)))); % 计算零交叉点
%             changePoints = zeroCrossings + 1; % 调整索引
%             lowerBound = max(changePoints - delta, 1); % 确保不小于0
%             upperBound = min(changePoints + delta, size(H1,3)); % 确保不大于endValue
%             allChangePoints = [allChangePoints; (lowerBound:upperBound)']; 
%         end
%     end
%     H2=imag(diff(H1));
%     for i =1:size(H2,1)
%         for j =1:size(H2,2)
%             zeroCrossings = find(diff(sign(H2(i,j,:)))); % 计算零交叉点
%             changePoints = zeroCrossings + 1; % 调整索引
%             lowerBound = max(changePoints - delta, 1); % 确保不小于0
%             upperBound = min(changePoints + delta, size(H1,3)); % 确保不大于endValue
%             allChangePoints = [allChangePoints;(lowerBound:upperBound)']; 
%         end
%     end
    
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
    
    
%     size(H2,3)
%     size(H1,3)
   
    t0 = clock;
    if(strcmp(func2str(func),'VF'))
        [G,W] = func(G,W,F2,H2,opts,valleypeak);
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
    
     
    [r2,~] = passivity_violation(G);
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

