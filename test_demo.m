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
opts.vf_niter2 =30;


opts.weightparam = 3;


opts.Sample= 1;%default%不需要采样设为2%%%%%基本没问题,三个参数delta delta1 delta2调一调效果就很好
% opts.sample_add = 0;
opts.enableSVD = 0;
opts.jk=1.4;%default%1.4
opts.final_test = 0;

windowSize = 30;
proximityThreshold = 20;
numRuns_start = 5;
numRuns_end = 5;

filenames = {'../../data/channel.S2P', ...
    '../../data/pll_sa_doubler_spur_ind.S5P', ...
    '../../data/pll_cko_lo1_top_l.s34p', ...
    '../../data/sp125_uniform.S64P', ...
    '../../data/pll_testcase.s138p', ...
    '../../data/graphs_typical_-40_0.s229p', ...
    '../../data/1848-191031.s384p'};

%% 是否要用到并行运算
% % 获取计算机的最大可用核心数（线程数）
% max_workers = feature('numCores'); % 返回可用核心数（包含超线程的数量）
% % 检查并行池是否已开启
% p = gcp('nocreate'); % 如果池未开启，则返回空
% if isempty(p)
%     % 开启并行池，设置 worker 数量为 max_workers
%     parpool('local', max_workers); % 最大并行 worker 数量
% end

% diary('eda_diary.txt');
% diary on;
for i =numRuns_start:length(filenames)
    fprintf("\n");
    t1=clock;
%     pool = gcp('nocreate'); % 获取当前并行池，如果没有则返回空
%     if ~isempty(pool)
%         delete(pool); % 删除并行池
%     end
    filename = filenames{i};  % 获取当前文件名
    [G,W,F,H,info]=pmm_S(filename,opts,windowSize);
    [residues,poles]=ss2pr(G.A,G.B,G.C);
    Hinf=G.D;

%% Enforce DC

   H_delta = G.D - G.C * (G.A \ G.B)  ;
   dc_error = max(vec(abs(H_delta - H(:,:,1))));
   if dc_error>1e-10
   all_error2=calculate(H,F,Hinf,residues,poles)
   
%    for ir = 1:length(poles)
%        H_delta = H_delta + residues(:,:,ir)./(- poles(ir));
%    end
   H_delta = real(H(:,:,1)-H_delta);
   alpha=-0.001;
   res   = -H_delta.*alpha; 
   poles = [alpha;poles];
   residues(:,:,2:size(residues,3)+1)=residues;
   residues(:,:,1)=res;
   end

%    e2-e1
    
   
   t2=clock;
   fprintf("非轻量化总程序耗时：%.5f s\n",etime(t2,t1));
   %% 计算误差
   all_error=calculate(H,F,Hinf,residues,poles)
%    delta=all_error2-all_error  
    
%     generate_model_dat(poles,residues,Hinf);
    
    if i == numRuns_end
        break;
    end
end
% diary off;
% 
% profile off
% 
% % 4. 生成包含时间戳的文件名
% filename = filenames{i};  % 输入文件名，可以修改为所需名称
% timestamp = datestr(now, 'yyyymmdd_HHMMSS');  % 获取当前时间，格式为 "yyyymmdd_HHMMSS"
% fullFilename = sprintf('%s_%s.mat', filename, timestamp);  % 组合文件名和时间戳
% 
% % 5. 获取性能分析数据并保存到文件中
% p = profile('info');  % 获取分析数据
% profile viewer
% save(fullFilename, 'p');  % 保存数据到 .mat 文件中
% 
% % 输出提示信息
% fprintf('性能分析数据已保存到文件：%s\n', fullFilename);
function all_error=calculate(H,F,Hinf,residues,poles)
    length(poles)
%          H2=[];
    for ik = 1:size(H,3)
        H2(:,:,ik) = Hinf;
        for ir = 1:length(poles)
        H2(:,:,ik) = H2(:,:,ik) + residues(:,:,ir)./(1j*2*pi*F(ik) - poles(ir));
        end
    end
    err1=0;
    err2=0;
    for w=1:size(H,3)
        err1 = err1 + norm(H(:,:,w) - H2(:,:,w),2);
        err2 = err2 + norm(H(:,:,w),2);
    end
%     for q=1:size(H,1)
%          for m=1:size(H,2)
%      err1 = err1 + norm(reshape(H(q,m,2:end) - H2(q,m,2:end), [], 1), 2);
%      err2 = err2 + norm(reshape(H(q,m,2:end), [], 1), 2);
%          end
%     end
    all_error=err1/err2*100;
    DC_ERROR =max(vec(abs(H2(:,:,1)-H(:,:,1))))
    
    %%Plot kind1
%     figure;
%     for idx = 1:4
%     subplot(2, 2, idx);
%     hold on;
%     % H1 和 H2 的实部和虚部曲线
%     plot(F, abs(reshape(H(idx,idx,:), [size(H,3), 1])), '-', 'Color', '#1f77b4', 'LineWidth', 1.5);
%     plot(F, abs(reshape(H2(idx,idx,:), [size(H,3), 1])),'-', 'Color', '#d62728', 'LineWidth', 1.5);
%     plot(F, abs(reshape(H(idx,idx,:) - H2(idx,idx,:), [size(H,3), 1]))*100,'-', 'Color', '#2ca02c', 'LineWidth', 1.5); 
%     title(['Subplot for H(', num2str(idx), ',', num2str(idx), ',:)']);
%     xlabel('Frequency(Hz)');
%     ylabel('Amplitude');
%     legend('H1', 'H2','Error(*100)');
%     % 设置期刊风格的颜色、字体
%     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
%     grid on;
%     hold off;
%     end
      
%      %%Plot kind2
%      % 计算 H 和 H2 的差的绝对值
%     abs_diff = abs(H - H2);
%     % 将差值拉伸为一维数组以供绘制
%     abs_diff_flat = abs_diff(:);
%     % 创建直方图
%     figure;
%     histogram(abs_diff_flat, 30, 'FaceColor', '#1f77b4', 'EdgeColor', '#d62728');
%     % 添加标题和标签
%     title('Error', 'FontSize', 14, 'FontName', 'Times New Roman');
%     xlabel('Magnitude', 'FontSize', 12, 'FontName', 'Times New Roman');
%     ylabel('Num', 'FontSize', 12, 'FontName', 'Times New Roman');
%     % 设置期刊风格
%     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
%     grid on;

      
     
end






