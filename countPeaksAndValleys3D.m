function [totalPeakNum] = countPeaksAndValleys3D(data, windowSize,proximityThreshold)
    % countPeaksAndValleys3D 统计三维复数矩阵中的波数量
    %
    % 用法:
    %   totalPeakNum = countPeaksAndValleys3D(data, windowSize)
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
                lastIndex = k;
            elseif (k - lastIndex) <= proximityThreshold
                mergedIndices(k) = false;
            else
                mergedIndices(k) = true;
                lastIndex = k;
            end
        end
    end
    
    % 计算总波峰和波谷数量
    totalPeakNum = round(sum(mergedIndices)/2);

    fprintf("time use %.4f seconds\n", toc);
end

