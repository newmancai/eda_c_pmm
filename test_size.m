A = [1, 2; 3, 0; 4, 5];  % Non-square matrix
block_size = 2;  % Size of identity matrix B

C_sparse = kron_optimized(A, block_size);  % Call the function
full(C_sparse)  % Convert the sparse matrix to full to view the result

disp(kron(A,eye(2)));
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

