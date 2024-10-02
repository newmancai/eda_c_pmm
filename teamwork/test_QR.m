
function test_QR()
    clc
    clear
    % 设置矩阵大小
    m = 127; % 行数
    n = 30; % 列数
    A = rand(m, n); % 生成随机矩阵 A

    % 测试 MATLAB 自带的 QR 分解
    tic; % 开始计时
    R = qr(A); % 计算经济型 QR 分解
    R = R(1:n, :); % 取前 n 行的 R
    R = triu(R);
    matlab_time = toc; % 结束计时
    fprintf('MATLAB QR time: %.6f seconds\n', matlab_time);

    tic; % 开始计时
    [Q,R11] = qr(A,0); % 计算经济型 QR 分解
    R11;
    matlab_time = toc; % 结束计时
    fprintf('MATLAB QR time: %.6f seconds\n', matlab_time);

end

function R = fastQR(A)
    % 快速 QR 分解的实现
    [m, n] = size(A);
    R = A;

    for k = 1:min(m-1, n)
        % 计算 Householder 反射
        x = R(k:m, k);
        e = zeros(size(x));
        e(1) = norm(x);  % 反射向量
        u = x - e;  % Householder 向量
        u = u / norm(u);  % 单位化
        
        % 更新矩阵 R
        R(k:m, k:n) = R(k:m, k:n) - 2 * (u * (u' * R(k:m, k:n)));
    end
end