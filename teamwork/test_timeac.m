close all
clear
clc
addpath('..');
pmm_setup

opts=pmm_default_S;
opts.method='vf_only';
opts.enforceDC = 1;
%opts.vf_niter1 =30;
opts.vf_niter2 =100;

windowSize = 20;
proximityThreshold =30;
numRuns = 1;

filenames = {'channel.s2p', 'pll_sa_doubler_spur_ind.s5p', 'sp125_uniform.s64p'};
for fileIndex = 1:length(filenames)
    filename = filenames{fileIndex};  % 获取当前文件名
    % 运行代码
    [G,W,F,H,info] = pmm_S(filename, opts, windowSize, proximityThreshold);
    print_info(info,extract_filename(filename));
    % 你可以在这里处理运行结果，比如存储结果或打印输出
end

%[G,W,F,H,info]=pmm_S(filename,opts,windowSize,proximityThreshold);
%[residues,poles]=ss2pr(G.A,G.B,G.C);
%Hinf=G.D;
%generate_model_dat(poles,residues,Hinf);
%figure; plot_xf(G,F,H);legend('Data','Model','Error');
%figure; plot_haeig(G);


function new_filename = extract_filename(original_filename)
    % 提取文件名前部分并保留后缀
    % 输入:
    %   original_filename - 原始文件名 (例如 'sp125_uniform.s64p')
    % 输出:
    %   new_filename - 提取后的文件名 (例如 'sp125.s64p')

    % 找到下划线的位置
    underscore_idx = strfind(original_filename, '_');

    % 使用 fileparts 分离文件名和文件后缀
    [~, name, ext] = fileparts(original_filename);

    % 如果找到下划线，提取下划线前面的部分，并加上文件后缀
    if ~isempty(underscore_idx)
        new_filename = [name(1:underscore_idx(1)-1), ext];
    else
        new_filename = [name, ext]; % 没有下划线的情况，直接使用原名和后缀
    end
end

function print_info(info,filename)
%sss
    col_widths = [15, 15, 15, 15, 15, 10]; % 各列宽度
    total_length = sum(col_widths) + length(col_widths) - 1; % 计算总长度（包含空格）
    fprintf('%-15s %-15s %-15s %-15s %-15s %-15s %-10s\n', ...
        'DocName','FuncName', 'Time', 'Error(Norm)', 'Error(DC)', 'K_accuracy', 'Passivity');
    fprintf('%s\n', repmat('-', 1, total_length));
    
    % 输出每个信息的内容
    for c = 1:length(info)
        errorWithPercent = sprintf('%.5f%%', info{c}.error * 100);
        fprintf('%-15s %-15s %-15.5f %-15s %-15.5f %-15.5f %-10s\n', ...
            filename,info{c}.func, info{c}.time, errorWithPercent, info{c}.dc_error, info{c}.k_accuracy, info{c}.passivity);
    end
end