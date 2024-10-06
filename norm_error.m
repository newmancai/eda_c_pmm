function [err]=norm_error(S,F,H1,H2)
    if nargin ==4
        if size(H2,1) == 1
            H2 = reshape(H2,[1,1,size(H2,2)]);
    elseif size(H1,1)~=size(H2,1)
            [n,~,m] = size(H1);       % 对称矩阵的大小
            H_new = zeros(n, n, m);   % 初始化nxnxm的矩阵

            for i = 1:m
                % 提取H2的第i列并将其重新构造为对称矩阵
                upper_tri = H2(:, i);   % 第 i 列包含对称矩阵的上三角部分
                symmetric_matrix = zeros(n, n);  % 用于存储对称矩阵

                % 填充上三角部分
                idx = 1;
                for row = 1:n
                    for col = row:n
                        symmetric_matrix(row, col) = upper_tri(idx);
                        symmetric_matrix(col, row) = upper_tri(idx); % 填充下三角部分
                        idx = idx + 1;
                    end
                end

                % 将对称矩阵放入 H_new 的第 i 层
                H_new(:, :, i) = symmetric_matrix;
            end
            H2 = H_new;
        end
    end
    if nargin < 4
        [residues,poles]=ss2pr(S.A,S.B,S.C);
        Hinf=S.D;
        M = size(residues,1);
        Nf = size(F,1); 
        [Nq, ~] = size(poles);  
        s = 2*pi*1i*F;
        H2 = zeros(M,M,Nf); 
        H2 = H2 + Hinf; 
        for k = 1:Nq
            factor = 1 ./ (s - poles(k));  % 201 x 1 矩阵
            H2(:,:, :) = H2(:,:, :) + residues(:,:,k) .* reshape(factor, 1, 1, Nf);
        end
    end
    %H2=ss_xf(S,F);
    if ndims(H2) == 3
        N = size(H2, 3);  % 获取频率点数量
        H_diff = H1 - H2;
        err_numerator = sum(arrayfun(@(p) norm(H_diff(:,:,p), 2), 1:N));
        err_denominator = sum(arrayfun(@(p) norm(H1(:,:,p), 2), 1:N));
        
    elseif ndims(H2) == 2
        N = size(H2, 2);  % 获取频率点数量
        err_numerator = 0;
        err_denominator = 0;
        for p = 1:N
            % 计算误差的分子部分: norm(F_p - S_p)
            err_numerator = err_numerator + norm(H1(:,p) - H2(:,p), 2);
            % 计算误差的分母部分: norm(S_p)
            err_denominator = err_denominator + norm(H1(:,p), 2);
        end
    end
    err = err_numerator / err_denominator;
end
