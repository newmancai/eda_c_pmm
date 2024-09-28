clc
clear
opts = pmm_default_S;
opts.enforceDC = 1;
fileList = {'HHM1506.s3p','IEEE39BUS1.s1p','channel.s2p'};  % 这里列出所有要读取的文件
n = length(fileList);  % 获取文件数量
for i = 1:n
    filename = fileList{i};  % 获取当前文件名
    [freq,Scatter_params,param_type] = readTouchstone(filename);
    [Frequencies,Parameters]=ldstone(filename,opts);
    if((compareComplexArrays(Frequencies,freq))&&(compareComplexArrays(Parameters,Scatter_params)))
         fprintf('读取%s参数的%s成功 \n', param_type,filename);
    else
         error('读取文件失败');
    end
end

function result = compareComplexArrays(A, B, tolerance)
    % compareComplexArrays 判断两个复数数组是否有 99% 或更多的元素在容差范围内相等
    %
    % 输入:
    %   A - 第一个复数数组
    %   B - 第二个复数数组
    %   tolerance - 容差值，用于判断相等的标准
    %
    % 输出:
    %   result - 布尔值，若 99% 或更多的元素在容差范围内相等则返回 true，否则返回 false

    % 确保 A 和 B 的维度一致
    if nargin < 3
        tolerance = 1e-5;  % 可以根据需要调整容差大小
    end

    if numel(A) ~= numel(B)
        error('数组 A 和 B 的大小不一致');
    end

    A_flat = A(:);
    B_flat = B(:);

    % 计算每个元素的绝对误差
    diff = abs(A_flat - B_flat);

    % 判断在容差范围内的元素数量
    within_tolerance = sum(diff < tolerance);

    % 计算满足条件的元素百分比
    percent_within_tolerance = (within_tolerance / numel(A_flat)) * 100;

    % 判断是否有99%的元素在容差范围内
    result = percent_within_tolerance >= 99;
end