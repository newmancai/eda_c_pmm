bigH(:,:,1) = [1 2 3;
               4 5 6;
               7 8 9];
           
bigH(:,:,2) = [10 11 12;
               13 14 15;
               16 17 18];
Nc = 3;
f_for = zeros(Nc*Nc, size(bigH, 3)); % 预分配
tell = 0;

for col = 1:Nc
    for row = 1:Nc
        tell = tell + 1;
        f_for(tell,:) = reshape(bigH(row, col, :), [1, size(bigH, 3)]);
    end
end

f_reshape = reshape(bigH, Nc*Nc, size(bigH, 3));
disp(f_for');
disp(f_reshape');
