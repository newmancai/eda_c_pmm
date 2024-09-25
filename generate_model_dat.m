function generate_model_dat(poles, residues, Hinf,team_name)
    % poles: 复数极点矩阵，大小为(Nq, 1)，每个极点为复数
    % residues: MxMxNq 的三维复数数组，每个极点对应一个 MxM 的复数留数矩阵
    % team_name: 队伍的名称，用于生成文件名。
    
    if nargin < 4
        team_name = 'niu';
    end

    % 参数校验
    if size(poles, 1) ~= size(residues, 3)
        error('极点的数量必须与留数矩阵的数量一致');
    end
    
    % 文件名生成
    filename = sprintf('%s_model.dat', team_name);
    fileID = fopen(filename, 'w');

    % 检查文件是否成功打开
    if fileID == -1
        error('无法打开文件: %s。请检查路径和写入权限。', filename);
    end
    
    % 计算极点信息
    Nq = size(poles, 1);  % 总极点数
    complex_poles = poles(imag(poles) ~= 0);  % 复极点
    real_poles = poles(imag(poles) == 0);     % 实极点
    Nqc = length(complex_poles) / 2;  % 复共轭极点对数
    Nqr = length(real_poles);         % 实极点数

    % 记录复极点和实极点的index
    complex_pole_indices = find(imag(poles) ~= 0);  % 复极点索引
    real_pole_indices = find(imag(poles) == 0);     % 实极点索引

    % 输出极点信息
    fprintf(fileID, 'Poles: %d %d %d\n', Nq, Nqc, Nqr);
    
    % 输出复极点
    for i = 1:2:2*Nqc
        fprintf(fileID, '%.16e %.16e\n', real(complex_poles(i)), imag(complex_poles(i)));
        fprintf(fileID, '%.16e %.16e\n', real(complex_poles(i+1)), imag(complex_poles(i+1)));
    end
    
    % 输出实极点
    for i = 1:Nqr
        fprintf(fileID, '%.16e\n', real(real_poles(i)));
    end
    
    % 输出留数信息
    M = size(residues, 1);  % 端口数 (矩阵大小)
    fprintf(fileID, 'Residues: %d\n', M);
    
    for k = 1:length(complex_pole_indices)
        residue_matrix = residues(:,:,complex_pole_indices(k));
         % 复留数矩阵，输出实部和虚部 (共 2M² 个实数)
        
         for row = 1:M
                for col = 1:M
                    fprintf(fileID, '%.16e %.16e ', real(residue_matrix(row, col)), imag(residue_matrix(row, col)));
                end        
        end
        fprintf(fileID, '\n');
    end
    
    for k = 1:length(real_pole_indices)
        residue_matrix = residues(:,:,real_pole_indices(k));
        % 实留数矩阵，仅输出实部 (共 M² 个实数)
        for row = 1:M
             for col = 1:M
                 fprintf(fileID, '%.16e ', real(residue_matrix(row, col)));
             end        
        end
        fprintf(fileID, '\n');
    end

    % 输出常数矩阵 Hinf，只输出实部
    fprintf(fileID, 'Hinf: %d\n', M);
    
    for row = 1:M
        fprintf(fileID, '%.16e ', real(Hinf(row, :)));
        fprintf(fileID, '\n');
    end
    
    % 关闭文件
    fclose(fileID);
end
