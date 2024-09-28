function [G,W]=VF(G0,W0,F,H,opts)
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
poles=[];
s1=2*pi*1j*F;
Gc=VFdriver(H,s1,poles,opt1,F); 
%读取opts的parametertype,若无定义则为Y
Gc.parametertype=optget(opts,'parametertype','S');

%如果一开始是passivity，则不拟合到f=0的点
G=ss_real(Gc);
[r2] = passivity_violation(G);

%if isempty(r2) && opts.enforce==1
if isempty(r2) ||(size(H,1)^2*opt1.N >65536)
    opts.enforceDC = 0;
    %G.D = real(H(:,:,1) + G.C * (G.A \ G.B)); %强制赋值
else
    opts.enforceDC = 1;
end

G = optimizeSystem(G,F,H,opts);
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
        C_sparse = kron_optimized(-G.B' * inv(G.A)', m);

        % Construct equality constraints for the optimization problem
        Aeq = [C_sparse, eye(m^2)];
        H0 = H(:, :, 1);
        beq = vec(H0);

        % Remove the first elements from F and H matrices
        F = F(2:end);
        H = H(:, :, 2:end);

        % Preprocess the inputs
        [R, QG, del2] = preprocess(F, H, G, opts);

        % Set up the quadratic programming problem
        Q = R' * R;
        f = -QG' * R;
        options = optimset('LargeScale', 'off', 'TypicalX', [vec(G.C); vec(G.D)], 'Display', 'off');

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

function C_sparse = kron_optimized(A, block_size)
    % A: input matrix (not necessarily square)
    % block_size: the size of the identity matrix B (B = eye(block_size))
    % C_sparse: the optimized sparse Kronecker product result

    [n_rows, n_cols] = size(A);  % Get the size of matrix A

    % Initialize arrays to store the row indices, column indices, and values
    row_idx = [];
    col_idx = [];
    val = [];
    
    % Fill the sparse Kronecker product matrix
    for i = 1:n_rows
        for j = 1:n_cols
            if A(i, j) ~= 0
                % Define the block indices corresponding to the Kronecker product
                row_start = (i - 1) * block_size + 1;
                row_end = i * block_size;
                col_start = (j - 1) * block_size + 1;
                col_end = j * block_size;
                
                % Append the block indices and values to the arrays
                row_idx = [row_idx; (row_start:row_end)'];
                col_idx = [col_idx; (col_start:col_end)'];
                val = [val; A(i, j) * ones(block_size, 1)];
            end
        end
    end
    
    % Construct the sparse matrix using the collected indices and values
    C_sparse = sparse(row_idx, col_idx, val, n_rows * block_size, n_cols * block_size);
end

