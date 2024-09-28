close all
addpath('..');
pmm_setup

opts=pmm_default_S;
opts.parametertype = 'Y';
opts.method='lc_dao';
[G,W,F,H]=pmm_S('pll_sa_doubler_spur_ind.S5P', 10, opts);
figure; plot_xf(G,F,H);legend('Data','Model','Error');
figure; plot_haeig(G);