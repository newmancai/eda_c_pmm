% % 检查并使用现有并行池，如果不存在则创建一个新的并行池
% tic;
if isempty(gcp('nocreate'))
    parpool(8); % 创建 8 个工作者的并行池
end
% parfor i = 1:8
%     %disp(['parfor Thread ', num2str(i), ' outputs: ', num2str(i)]);
% end
% parfor_time = toc;
% delete(gcp('nocreate'));
% disp(toc)
% 
% a = randn([2048,2048,4]);
% 
% tic
% for i = 1:4 % 串行计算
%     svd(a(:,:,i));
% end
% toc
% 
% tic
% parfor i = 1:4 % 并行计算
%    svd(a(:,:,i));
% end
% toc

% 初始化一个矩阵 A 和 B
n = 10;
N = 5;
result = zeros(n, N+1, N+1); % 预分配一个三维矩阵

% 使用 parfor 循环对每个 n 的切片进行赋值
parfor i = 1:n
    % 对 result(i,:,:) 进行操作，例如填充随机值
    result(i, :, :) = rand(N+1, N+1); % 随机生成一个 (N+1) x (N+1) 矩阵并赋值给切片
end

% 输出结果
disp(result);