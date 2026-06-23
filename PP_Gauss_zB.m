% PHUONG PHAP GAUSS GIAI PHUONG TRINH MA TRAN AX = B (B co z cot)
% ==========================================
% B1: Cau hinh du lieu
% ==========================================
A = [ 1  2  0  3 ;
      2  4  0  6 ;
      0  0  0  1 ;
      0  0  0  2 ;
      0  0  0  3 ];

B = [ 1   2 ;
      2   4 ;
      5   1 ;
      10  2 ;
      15  3 ];

An = [A, B]; % Ma trận bổ sung [A | B]

[m, cols] = size(An);
z = 2; % Số cột của ma trận B
n = cols - z; % Số ẩn (Số cột của A)

fprintf('Ma tran ban dau An (%d phuong trinh, %d an, B co %d cot):\n', m, n, z);
disp(An);

i = 1; j = 1; buoc = 1;

% ==========================================
% B2: Qua trinh thuan (Khu Gauss voi Partial Pivoting)
% ==========================================
while (i <= m) && (j <= n)
    % Tim phan tu co gia tri tuyet doi lon nhat trong cot j tu hang i den m
    [val, relative_p] = max(abs(An(i:m, j)));
    p = relative_p + i - 1;

    % Neu cot j toan so 0 (hoac gan bang 0), nhay sang cot tiep theo
    if val < 1e-10
        fprintf('\n[THONG BAO] Cot %d tai hang %d khong co Pivot, nhay sang cot %d...\n', j, i, j+1);
        j = j + 1;
        continue;
    end

    % Doi hang i voi hang p (hang chua phan tu troi)
    if p ~= i
        An([i, p], :) = An([p, i], :);
        fprintf('\n--- Da doi hang %d voi hang %d ---\n', i, p);
    end

    % --- IN VI TRI KHU VA PIVOT ---
    fprintf('\n[BUOC %d] -------------------------------------------\n', buoc);
    fprintf('Vi tri Pivot dang xet: (hang %d, cot %d)\n', i, j);
    fprintf('Gia tri Pivot: %.4f\n', An(i, j));

    % Tien hanh khu cac hang duoi hang i
    for l = i + 1 : m
        heso = An(l, j) / An(i, j);
        An(l, :) = An(l, :) - heso * An(i, :);
    end

    % --- IN MA TRAN SAU KHI KHU TAI BUOC NAY ---
    fprintf('Ma tran sau khi khu tai cot %d:\n', j);
    disp(An);

    i = i + 1; j = j + 1;
    buoc = buoc + 1;
end

fprintf('============================================\n');
disp('MA TRAN KET QUA (DANG BAC THANG):');
disp(An);

% ==========================================
% PHAN BIEN LUAN (Cho phuong trinh ma tran)
% ==========================================

id = zeros(1, m);
for r = 1:m
    idx = find(abs(An(r, :)) > 1e-10, 1);
    if ~isempty(idx)
        id(r) = idx;
    end
end

i_max = find(id ~= 0, 1, 'last');
if isempty(i_max); i_max = 0; end

if i_max == 0
    fprintf('B6: He co ma tran nghiem X tuy y tren khong gian R^(%dx%d)\n', n, z);
    return;
end

fprintf('B7: Bien luan nghiem (imax = %d, id(imax) = %d, n = %d):\n', i_max, id(i_max), n);

if id(i_max) > n
    fprintf('==> KET LUAN: PHUONG TRINH MA TRAN VO NGHIEM.\n');
else
    k = n - i_max;
    all_vars = 1:n;
    main_vars = id(1:i_max);
    pa = setdiff(all_vars, main_vars);

    if k == 0
        fprintf('==> KET LUAN: CO MA TRAN NGHIEM DUY NHAT.\n');
        X = zeros(n, z);
        for row = i_max:-1:1
            col_pivot = id(row);
            X(col_pivot, :) = (An(row, n+1:end) - An(row, col_pivot+1:n) * X(col_pivot+1:n, :)) / An(row, col_pivot);
        end
        disp('Ma tran nghiem X ='); disp(X);
    else
        fprintf('==> KET LUAN: CO VO SO MA TRAN NGHIEM.\n');
        fprintf('(Moi cot cua X se phu thuoc %d tham so).\n', k);
        fprintf('Cac hang cua ma tran X la an tu do: '); disp(pa);
    end
end
