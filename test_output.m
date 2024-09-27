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
original_filename1 = 'sp125_uniform.s64p';
original_filename2 = 'pll_sa_doubler_spur_ind.s5p';

new_filename1 = extract_filename(original_filename1);
new_filename2 = extract_filename(original_filename2);

disp(new_filename1); % 输出: sp125.s64p
disp(new_filename2); % 输出: pll.s5p

function new_filename = extract_filename(original_filename)
    % 提取文件名前部分并保留后缀
    % 输入:
    %   original_filename - 原始文件名 (例如 'sp125_uniform.s64p')
    % 输出:
    %   new_filename - 提取后的文件名 (例如 'sp125.s64p')

    % 找到下划线的位置
    underscore_idx = strfind(original_filename, '_');

    % 使用 fileparts 分离文件名和文件后缀
    [~, name, ext] = fileparts(original_filename);

    % 如果找到下划线，提取下划线前面的部分，并加上文件后缀
    if ~isempty(underscore_idx)
        new_filename = [name(1:underscore_idx(1)-1), ext];
    else
        new_filename = [name, ext]; % 没有下划线的情况，直接使用原名和后缀
    end
end