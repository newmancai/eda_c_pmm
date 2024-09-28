function [ valleypeakIndices,totalPeakNum] = countPeaksAndValleys3D(data, windowSize,proximityThreshold)
% countPeaksAndValleys3D 统计三维复数矩阵中的波数量
%
% 用法:
%   totalPeakNum = countPeaksAndValleys3D(data, windowSize, proximityThreshold)
%
% 输入参数:
%   data        - 一个 M*M*n 的三维复数矩阵（对称阵）
%   windowSize  - 用于计算局部均值和标准差的窗口大小
%
% 输出参数:
%   totalPeakNum - 总的波数量（经过调整）
tic
assert(ndims(data) == 3, '输入数据必须是一个三维矩阵');

[M, ~, n] = size(data);
absData = abs(data);

allExtremaIndices = false(n, 1);
tempExtremaIndices = [];
valleypeakIndices = [];

for idx = 1:(M * (M + 1) / 2)
    [i, j] = ind2sub([M, M], idx);
    slice = squeeze(absData(i, j, :));

    localMean = movmean(slice, windowSize, 'Endpoints', 'shrink');
    localStd = movstd(slice, windowSize, 'Endpoints', 'shrink');

    isPeak = (slice(2:n-1) > slice(1:n-2)) & (slice(2:n-1) > slice(3:n)) & ...
        (slice(2:n-1) > localMean(2:n-1) + localStd(2:n-1));
    isValley = (slice(2:n-1) < slice(1:n-2)) & (slice(2:n-1) < slice(3:n)) & ...
        (slice(2:n-1) < localMean(2:n-1) - localStd(2:n-1));
    allExtremaIndices(2:n-1) = allExtremaIndices(2:n-1) | isPeak | isValley;
end

% 合并相邻且接近的波峰和波谷
mergedIndices = false(n, 1);
lastIndex = 0;
for k = 1:n
    if allExtremaIndices(k)
        if lastIndex == 0
            mergedIndices(k) = true;
            valleypeakIndices=[valleypeakIndices;k];
            lastIndex = k;
        elseif (k - lastIndex) <= proximityThreshold
            mergedIndices(k) = false;
        else
            mergedIndices(k) = true;
             valleypeakIndices=[valleypeakIndices;k];
             lastIndex = k;
        end
    end
end

% 计算总波峰和波谷数量
totalPeakNum = round(sum(mergedIndices)/2)+1;

fprintf("time use %.4f seconds\n", toc);
% fprintf("第一个结果 %d ", totalPeakNum);


% new_data = ifft(data(valleypeakIndices));
% % new_data = [real(new_data);imag(new_data)];
% 
% [A, B] = prony(data(valleypeakIndices), M, M) % A 和 B 是 Prony 方法的输出

% % 进行 SVD
% [U, S, V] = svd(new_data);
% % 打印结果
% disp('U矩阵：');
% disp(U);
% disp('奇异值（对角矩阵形式）：');
% disp(S);
% disp('V矩阵：');
% disp(V);

end

