% 示例：分析 myScript.m 的性能并将结果存档
clc
clear
% 1. 启动性能分析
addpath('..');
profile on

% 2. 执行你想要分析的 M 文件或函数
[freq,Scatter_params,param_type,impedance] = readTouchstone('../../data/pll_sa_doubler_spur_ind.S5P');  % 或者调用你想要分析的函数
%[freq,Scatter_params,param_type,impedance] = readTouchstone('../../data/1848-191031.s384p');  % 或者调用你想要分析的函数

% 3. 停止性能分析
profile off

% 4. 生成包含时间戳的文件名
filename = 'profile_data';  % 输入文件名，可以修改为所需名称
timestamp = datestr(now, 'yyyymmdd_HHMMSS');  % 获取当前时间，格式为 "yyyymmdd_HHMMSS"
fullFilename = sprintf('%s_%s.mat', filename, timestamp);  % 组合文件名和时间戳

% 5. 获取性能分析数据并保存到文件中
p = profile('info');  % 获取分析数据
profile viewer
save(fullFilename, 'p');  % 保存数据到 .mat 文件中

% 输出提示信息
fprintf('性能分析数据已保存到文件：%s\n', fullFilename);
