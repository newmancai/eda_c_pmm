% 参数初始化
port_num = 3;  % 端口数量
n = 4;  % G.C 的行数
m = 5;  % G.C 的列数

% 随机生成 G.C 矩阵和 V_effective 矩阵
G.C = randn(n, m);  % G.C 为 n x m 矩阵
V_effective = randn(port_num * port_num, n);  % V_effective 为 P^2 x n 矩阵
I_P = eye(port_num);  % P x P 单位矩阵

%% 方法1：显式计算 Kronecker 积
% 计算 Psi
Psi = [];
for j = 1:port_num
    V_j = V_effective((j-1)*port_num + 1:j*port_num, :);
    Psi = [Psi, V_j];
end

% 计算 Kronecker 积
C1 = Psi * kron(I_P, G.C);

%% 方法2：逐块乘法计算，避免显式构造 Kronecker 积
C2 = zeros(size(Psi, 1), size(G.C, 2) * port_num);  % 预分配 C2 矩阵
for j = 1:port_num
    V_j = V_effective((j-1)*port_num + 1:j*port_num, :);  % 提取 P 行
    C2 = C2 + V_j * G.C;  % 块乘法并累加
end

%% 对比结果
disp('方法1 (显式 Kronecker 积) 计算的 C 矩阵：');
disp(C1);

disp('方法2 (逐块乘法) 计算的 C 矩阵：');
disp(C2);

% 检查两个结果是否相等
if norm(C1 - C2, 'fro') < 1e-10
    disp('两种方法的结果一致！');
else
    disp('两种方法的结果不一致。');
end