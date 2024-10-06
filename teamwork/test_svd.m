% 生成一个随机对称矩阵 A
n = 5; % 矩阵的维度
A = randn(n); 
A = (A + A') / 2; % 确保 A 是对称矩阵

% 对称矩阵 A 的 SVD 分解
S_A = svds(A,2); % S_A 中的对角线是奇异值
singular_values_A = diag(S_A); % 提取奇异值

% 取对称矩阵 A 的上三角部分
A_upper = triu(A);

% 上三角矩阵的 SVD 分解
S_upper = svds(A_upper,2); 
singular_values_upper = diag(S_upper); % 提取奇异值

% 比较奇异值（特征值的绝对值）
disp('对称矩阵 A 的奇异值（特征值的绝对值）：');
disp(singular_values_A);

disp('上三角矩阵的奇异值：');
disp(singular_values_upper);

% 结果判断
if isequal(round(singular_values_A, 10), round(singular_values_upper, 10))
    disp('奇异值相等！');
else
    disp('奇异值不相等。');
end
