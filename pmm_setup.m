function setup_pmm()
    % 测试 setup_path 的时间
    tic; % 开始计时
    setup_path; % 设置路径
    time_setup_path = toc; % 结束计时并获取时间
    fprintf('setup_path 用时: %.4f 秒\n', time_setup_path);
    
    % 测试 setup_cvx 的时间
    tic; % 开始计时
    %setup_cvx; % 设置 CVX
    time_setup_cvx = toc; % 结束计时并获取时间
    fprintf('setup_cvx 用时: %.4f 秒\n', time_setup_cvx);
    
    % 测试 setup_pso 的时间（如果需要的话）
    % tic; % 开始计时
    % setup_pso; % 设置 PSO
    % time_setup_pso = toc; % 结束计时并获取时间
    % fprintf('setup_pso 用时: %.4f 秒\n', time_setup_pso);
    
    % 测试 preinstall_methods 的时间
    tic; % 开始计时
    preinstall_methods; % 预安装方法
    time_preinstall_methods = toc; % 结束计时并获取时间
    fprintf('preinstall_methods 用时: %.4f 秒\n', time_preinstall_methods);

function preinstall_methods()
pmm_install('VF', 'Vector Fitting', 0, 0);
pmm_install('SDP', 'SDP Method', 0, 1);
pmm_install('ASYM', 'Asymptotic passivity', 0, 1);
pmm_install('LC', 'Local Compensation', 0, 1);
pmm_install('EPM', 'Eigenvalue Perturbation', 0, 1);
pmm_install('FRP', 'FRP Method', 0, 1);
pmm_install('DAO', 'Domain Alternated Optimization', 1, 0);
pmm_install('EPM2', 'Eigenvalue Perturbation improved', 0, 1);
pmm_install('EPM3', 'Eigenvalue Perturbation improved', 0, 1);
pmm_install('EPM4', 'Eigenvalue Perturbation improved', 0, 1);

function setup_path()
[root,name,ext]=fileparts(mfilename('fullpath'));
addpath(root);
addpath(sprintf('%s/utils',root));
addpath(sprintf('%s/methods',root));
addpath(sprintf('%s/methods/dao',root));
tic
addpath(sprintf('%s/methods/epm',root));
toc
tic
addpath(sprintf('%s/others/mfit',root));
toc
tic
addpath(sprintf('%s/others/cvx',root));
toc

function setup_pso()
[root,name,ext]=fileparts(mfilename('fullpath'));
dir=pwd;
cd(sprintf('%s/others/PSOt',root));
pso_setup;
cd(dir);
cvx_loaded=1;


function setup_cvx()
   
[root,name,ext]=fileparts(mfilename('fullpath'));
dir=pwd;
cd(sprintf('%s/others/cvx',root));
cvx_setup;
cd(dir);
cvx_loaded=1;
