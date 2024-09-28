function pmm_setup()
    setup_path; % 设置路径
    
    %setup_cvx; % 设置 CVX
    
    % setup_pso; % 设置 PSO

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
[root]=fileparts(mfilename('fullpath'));
addpath(root);
addpath(sprintf('%s/utils',root));
addpath(sprintf('%s/methods',root));
%addpath(sprintf('%s/methods/dao',root));
addpath(sprintf('%s/methods/epm',root));
addpath(sprintf('%s/others/mfit',root));
%addpath(sprintf('%s/others/cvx',root));

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
