function [err]=norm_error(S,F,H1)
    [residues,poles]=ss2pr(S.A,S.B,S.C);
    Hinf=S.D;
    M = size(residues,1);
    Nf = size(F,1); 
    [Nq, ~] = size(poles);  
    s = 2*pi*1i*F;
    H2 = zeros(M,M,Nf); 
    H2 = H2 + Hinf; 
    % 第二步：加实极点的贡献
    for k = 1:Nq
        for i = 1:Nf
            H2(:,:,i) = H2(:,:,i) + residues(:,:,k) ./ (s(i) - poles(k));
        end
    end
    %H2=ss_xf(S,F);
    N = size(H2, 3);  % 获取频率点数量
    err_numerator = 0;
    err_denominator = 0;
    for p = 1:N
        % 计算误差的分子部分: norm(F_p - S_p)
        err_numerator = err_numerator + norm(H1(:,:,p) - H2(:,:,p), 2);
        % 计算误差的分母部分: norm(S_p)
        err_denominator = err_denominator + norm(H1(:,:,p), 2);
    end
    err = err_numerator / err_denominator;
end
