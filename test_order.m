opts=pmm_default_S;
filename = 'channel.s2p';
%filename = 'pll_sa_doubler_spur_ind.s5p';
%filename = 'sp125_uniform.s64p';
[F,H] = readTouchstone(filename,opts);

windowSize = 20; 
proximityThreshold = 30;
[totalPeakNum] = countPeaksAndValleys3D(H, windowSize,30);

fprintf('总波数量: %d\n', totalPeakNum);

