function [G,W]=VF(G0,W0,F,H,opts,valleypeak)
% VF -- generate state-space model from tabulated data using 
% Vector fitting method.
% This is simply an interface to VFdriver from Matrix Fitting Toolbox,
% This toolbox can be downloaded from 
% http://www.energy.sintef.no/produkt/VECTFIT/index.asp
%
% Copytright Zuochang Ye, 2012

opt1.N=opts.q ;%           %Order of approximation. 
opt1.Niter1=optget(opts,'vf_niter1', 200);    %Number of iterations for fitting sum of elements (fast!) --> Improved initial poles
opt1.Niter2=optget(opts,'vf_niter2', 100);    %Number of iterations for matrix fitting 
opt1.asymp=2;      %Fitting includes D   
opt1.logx=0;       %=0 --> Plotting is done using linear abscissa axis 
opt1.poletype=optget(opts,'poletype','lincmplx');
opt1.spy1=0; 
opt1.spy2=0; 
opt1.logx=0; 
opt1.logy=1; 
opt1.errplot= 1;  %%%%%%%%%%%%
opt1.phaseplot=0;
opt1.screen=0; %optget(opts,'verbose',0);
opt1.plot=opts.plot;
opt1.cmplx_ss=0;
opt1.weightparam=optget(opts,'wgt_scheme',3); %weight(s)=1/sqrt(abs(Hij(s)));
opt1.enableSVD = optget(opts,'enableSVD',0);
poles=[];
s1=2*pi*1j*F;
if ~opt1.enableSVD
    Gc=VFdriver(H,s1,poles,opt1,F,valleypeak); 
else
    Gc=VFdriver_SVD(H',s1,poles,opt1,F,valleypeak); 
end
%读取opts的parametertype,若无定义则为Y
Gc.parametertype=optget(opts,'parametertype','S');

%如果一开始是passivity，则不拟合到f=0的点
G=ss_real(Gc);
% [r2] = passivity_violation(G);

%if isempty(r2) && opts.enforce==1
if size(H,1)^2*opt1.N >65536
    opts.enforceDC = 0;
    %G.D = real(H(:,:,1) + G.C * (G.A \ G.B)); %强制赋值
else
    opts.enforceDC = 1;
end

%G = optimizeSystem(G,F,H,opts);
W=[];

end

function G = optimizeSystem(G, F, H, opts)
    % Input:
    % G - System structure containing matrices A, B, C, D
    % F - Input matrix for preprocessing
    % H - Output matrix for preprocessing
    % opts - Option structure

    % Check if F(1) is 0 and enforceDC option is enabled
    if F(1) == 0 && optget(opts, 'enforceDC', 0)
        [n, m] = size(G.B);  % Get the dimensions of G.B

        % Generate sparse matrix using Kronecker product
        C_sparse = kron(-G.B' * inv(G.A)', eye(m));

        % Construct equality constraints for the optimization problem
        Aeq = [C_sparse, eye(m^2)];
        H0 = H(:, :, 1);
        beq = vec(H0);

        % Remove the first elements from F and H matrices
        F = F(2:end);
        H = H(:, :, 2:end);

        % Preprocess the inputs
        [R, QG, ~] = preprocess(F, H, G, opts);

        % Set up the quadratic programming problem
        Q = R' * R;
        f = -QG' * R;
        vecG = [vec(G.C); vec(G.D)];
        % samll num
        small_value = 1e-8;  % 可以根据需要调整这个值

        % replace 0
        vecG(vecG == 0) = small_value;
        options = optimset('LargeScale', 'off', 'TypicalX', vecG, 'Display', 'off');

        % Solve the quadratic programming problem
        y = quadprog(Q, f, [], [], Aeq, beq, [], [], [], options);

        % Take the real part of the result
        y = real(y);

        % Reshape the solution into G.C and G.D matrices
        G.C = reshape(y(1:m * n), m, n);
        G.D = reshape(y((m * n + 1):(n * m + m^2)), m, m);

        % Update G.D with the current values of G.C, G.A, and G.B
        G.D = H0 + G.C * (G.A \ G.B);
    end
end

