close all
clear
clc

if strcmp(profile('status').ProfilerStatus,'on')  
    profile off;  
end
profile on
disp('性能分析工具已开启。');

addpath('..');
pmm_setup

opts=pmm_default_S;
opts.method='vf_only';
opts.enforceDC = 1;
opts.vf_niter1 =300;
opts.vf_niter2 =10;


opts.weightparam = 3;


opts.Sample= 1;%default%不需要采样设为2%%%%%基本没问题,三个参数delta delta1 delta2调一调效果就很好
opts.jk=1.4;%default%1.4

windowSize = 30;
proximityThreshold = 20;
numRuns_start = 2;
numRuns_end = 2;

filenames = {'../../data/pll_sa_doubler_spur_ind.S5P','../../data/pll_cko_lo1_top_l.s34p','../../data/pll_testcase.s138p','../../data/sp125_uniform.S64P','../../data/1848-191031.s384p', '../../data/channel.S2P', '../../data/graphs_typical_-40_0.s229p'};

for i =numRuns_start:length(filenames)
    t1=clock;
    pool = gcp('nocreate'); % 获取当前并行池，如果没有则返回空
    if ~isempty(pool)
        delete(pool); % 删除并行池
    end
    filename = filenames{i};  % 获取当前文件名
    [G,W,F,H,info]=pmm_S(filename,opts,windowSize);
    [residues,poles]=ss2pr(G.A,G.B,G.C);
    Hinf=G.D;
    generate_model_dat(poles,residues,Hinf);
    t2=clock;
    fprintf("非轻量化总程序耗时：%.5f s\n",etime(t2,t1))
    if i == numRuns_end
        break;
    end
end

profile off

% 4. 生成包含时间戳的文件名
filename = filenames{i};  % 输入文件名，可以修改为所需名称
timestamp = datestr(now, 'yyyymmdd_HHMMSS');  % 获取当前时间，格式为 "yyyymmdd_HHMMSS"
fullFilename = sprintf('%s_%s.mat', filename, timestamp);  % 组合文件名和时间戳

% 5. 获取性能分析数据并保存到文件中
p = profile('info');  % 获取分析数据
profile viewer
save(fullFilename, 'p');  % 保存数据到 .mat 文件中

% 输出提示信息
fprintf('性能分析数据已保存到文件：%s\n', fullFilename);
