function [err]=norm_error(S,F,H1)
    H2=ss_xf(S,F);
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
