close all
clear
clc


addpath('..');
pmm_setup

opts=pmm_default_S;
opts.method='vf_only';
opts.enforceDC = 1;
opts.vf_niter1 =300;
opts.vf_niter2 =1000;


opts.weightparam = 3;


opts.Sample= 1;%default%不需要采样设为2%%%%%基本没问题,三个参数delta delta1 delta2调一调效果就很好
opts.jk=1.4;%default%1.4

windowSize = 20;
proximityThreshold = 30;
numRuns = 1;

filenames = {'pll_cko_lo1_top_l.s34p','pll_sa_doubler_spur_ind.S5P','pll_testcase.s138p','sp125_uniform.S64P','1848-191031.s384p', 'channel.S2P', 'graphs_typical_-40_0.s229p'};

for i =1:length(filenames)
    t1=clock;
    pool = gcp('nocreate'); % 获取当前并行池，如果没有则返回空
    if ~isempty(pool)
        delete(pool); % 删除并行池
    end
    filename = filenames{i};  % 获取当前文件名
    [G,W,F,H,info]=pmm_S(filename,opts,windowSize,proximityThreshold);
    print_info(info,extract_filename(filename));
    [residues,poles]=ss2pr(G.A,G.B,G.C);
    Hinf=G.D;
%     order = process_frequency_response(H,F*2*pi*1j, poles, residues,Hinf);
%     fprintf("推荐阶数：%d \n",order);
%     if (order*size(H,1)^2<7e4)&&(order>size(poles,1))
%         [G1,~,F1,H1,info1]=pmm_S(filename,opts,windowSize,proximityThreshold,order);
%         print_info(info1,extract_filename(filename));
%         [residues1,poles1]=ss2pr(G1.A,G1.B,G1.C);
%         Hinf1=G1.D;
%     end
    generate_model_dat(poles,residues,Hinf);
%     if ~isempty(info1)&&(info1{1}.k_accuracy>info{1}.k_accuracy)
%         generate_model_dat(poles1,residues1,Hinf1);
%     else
%         generate_model_dat(poles,residues,Hinf);
%     end
    t2=clock;
    fprintf("非轻量化总程序耗时：%.5f s\n",etime(t2,t1))
    % [residues,poles]=ss2pr(G.A,G.B,G.C);
    % Hinf=G.D;
    % generate_model_dat(poles,residues,Hinf);
    % figure; plot_xf(G,F,H);legend('Data','Model','Error');
    % figure; plot_haeig(G);
end


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
    col_widths = [15,15, 15, 15, 15, 15, 10]; % 各列宽度
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