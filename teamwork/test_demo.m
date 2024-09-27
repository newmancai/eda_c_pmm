close all
clear
addpath('..');
pmm_setup

opts=pmm_default_S;

%filename = 'channel.s2p';
%filename = 'pll_sa_doubler_spur_ind.s5p';
filename = 'sp125_uniform.s64p';

opts.method='vf_only';
opts.enforceDC = 1;
%opts.vf_niter1 =30;
%opts.vf_niter2 =5;
windowSize = 20;
proximityThreshold =30;
[G,W,F,H,info]=pmm_S(filename,opts,windowSize,proximityThreshold);
[residues,poles]=ss2pr(G.A,G.B,G.C);
Hinf=G.D;
generate_model_dat(poles,residues,Hinf);
figure; plot_xf(G,F,H);legend('Data','Model','Error');
figure; plot_haeig(G);