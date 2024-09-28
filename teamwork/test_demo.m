close all
clear
clc
<<<<<<< HEAD
t1=clock;
=======
>>>>>>> 9795d410fe9d13233d9941c3bede4a0575432838
addpath('..');
pmm_setup
opts=pmm_default_S;
opts.method='vf_only';
opts.enforceDC = 1;
<<<<<<< HEAD
opts.enforce = 1;
opts.Sample =1;
opts.vf_niter1 =600;
opts.vf_niter2 =300;

windowSize = 30;
=======
opts.vf_niter1 =300;
opts.vf_niter2 =1000;

windowSize = 20;
>>>>>>> 9795d410fe9d13233d9941c3bede4a0575432838
proximityThreshold =30;
numRuns = 1;

filenames = {'channel.s2p', 'pll_sa_doubler_spur_ind.s5p', 'sp125_uniform.s64p'};

filename = filenames{1};  % 获取当前文件名
<<<<<<< HEAD

[G,~,F,H,info]=pmm_S(filename,opts,windowSize,proximityThreshold);
print_info(info,extract_filename(filename));
[residues,poles]=ss2pr(G.A,G.B,G.C);
Hinf=G.D;
order = process_frequency_response(H,F*2*pi*1j, poles, residues,Hinf,0);
fprintf("推荐阶数：%d ",order);
if (order*size(H,1)^2<7e4)
    [G1,~,F1,H1,info1]=pmm_S(filename,opts,windowSize,proximityThreshold,order);
    print_info(info1,extract_filename(filename));
    [residues1,poles1]=ss2pr(G.A,G.B,G.C);
    Hinf1=G.D;
end
if(info1{1}.k_accuracy>info{1}.k_accuracy)
    generate_model_dat(poles,residues,Hinf);
else
    generate_model_dat(poles1,residues1,Hinf1);
end
%figure; plot_xf(G,F,H);legend('Data','Model','Error');
%figure; plot_haeig(G);
t2=clock;
etime(t2,t1)
=======
[G,W,F,H,info]=pmm_S(filename,opts,windowSize,proximityThreshold);
print_info(info,extract_filename(filename));
[residues,poles]=ss2pr(G.A,G.B,G.C);
Hinf=G.D;
generate_model_dat(poles,residues,Hinf);
figure; plot_xf(G,F,H);legend('Data','Model','Error');
figure; plot_haeig(G);

>>>>>>> 9795d410fe9d13233d9941c3bede4a0575432838

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
<<<<<<< HEAD
    col_widths = [15,15, 15, 20, 15, 15, 10]; % 各列宽度
    total_length = sum(col_widths) + length(col_widths) - 1; % 计算总长度（包含空格）
    fprintf('%-15s %-15s %-15s %-20s %-15s %-15s %-10s\n', ...
=======
    col_widths = [15, 15, 15, 15, 15, 10]; % 各列宽度
    total_length = sum(col_widths) + length(col_widths) - 1; % 计算总长度（包含空格）
    fprintf('%-15s %-15s %-15s %-15s %-15s %-15s %-10s\n', ...
>>>>>>> 9795d410fe9d13233d9941c3bede4a0575432838
        'DocName','FuncName', 'Time', 'Error(Norm)', 'Error(DC)', 'K_accuracy', 'Passivity');
    fprintf('%s\n', repmat('-', 1, total_length));
    
    % 输出每个信息的内容
    for c = 1:length(info)
<<<<<<< HEAD
        errorWithPercent = sprintf('%.9f%%', info{c}.error * 100);
        fprintf('%-15s %-15s %-15.5f %-20s %-15.5f %-15.8f %-10s\n', ...
=======
        errorWithPercent = sprintf('%.5f%%', info{c}.error * 100);
        fprintf('%-15s %-15s %-15.5f %-15s %-15.5f %-15.5f %-10s\n', ...
>>>>>>> 9795d410fe9d13233d9941c3bede4a0575432838
            filename,info{c}.func, info{c}.time, errorWithPercent, info{c}.dc_error, info{c}.k_accuracy, info{c}.passivity);
    end
end