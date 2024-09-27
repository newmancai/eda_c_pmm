% 示例数据
info = {
    struct('func', 'VF', 'time', 15.34375, 'error', 0.06340, 'dc_error', 3.088282e-18, 'k_accuracy', 0.1267977, 'passivity', 'non-passive');
    struct('func', 'AF', 'time', 12.45678, 'error', 0.04567, 'dc_error', 1.234567e-17, 'k_accuracy', 0.0987654, 'passivity', 'passive')
};

% 输出表头
fprintf('%-15s %-15s %-15s %-15s %-15s %-10s\n', ...
    'FuncName', 'Time', 'Error(Norm)', 'Error(DC)', 'K_accuracy', 'Passivity');
fprintf('------------------------------------------------------------------------------------------\n');

% 输出每一行数据
for c = 1:length(info)
    % 将error值和百分号结合成字符型
    errorWithPercent = sprintf('%.1f%%', info{c}.error * 100);
    
    fprintf('%-15s %-15.5f %-15s %-15.5f %-15.5f %-10s\n', ...
        info{c}.func, info{c}.time, errorWithPercent, info{c}.dc_error, info{c}.k_accuracy, info{c}.passivity);
end