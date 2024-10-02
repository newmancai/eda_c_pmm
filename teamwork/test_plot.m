opts=pmm_default_S;
filename = '../../data/pll_sa_doubler_spur_ind.s5p';
%filename = 'channel.s2p';
[F,H] = readTouchstone(filename,opts);
% 获取矩阵尺寸
[M, ~, n] = size(H);

% 对 H 取以10为底的对数
logH = log10(H);

% 创建图形窗口
figure;

% 绘制 H 的实部（log10后）
for i = 1:M
    for j = 1:M
        subplot(M, M, (i-1)*M + j);  % 创建子图
        plot(F, abs(squeeze(H(i,j,:))));
        title(sprintf('Log10(H(%d,%d)) abs Part', i, j));
        xlabel('Frequency');
        ylabel('abs Part');
    end
end

% % 新建图形窗口以绘制虚部（log10后）
% figure;
% 
% % 绘制 H 的虚部（log10后）
% for i = 1:M
%     for j = 1:M
%         subplot(M, M, (i-1)*M + j);  % 创建子图
%         plot(F, imag(squeeze(logH(i,j,:))));
%         title(sprintf('Log10(H(%d,%d)) Imaginary Part', i, j));
%         xlabel('Frequency');
%         ylabel('Imaginary Part');
%     end
% end
