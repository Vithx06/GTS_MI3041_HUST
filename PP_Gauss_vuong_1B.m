% PHUONG PHAP GAUSS GIAI HE PHUONG TRINH (CO BIEN LUAN NGHIEM)
% B1: Cau hinh du lieu
An = [
 1  2  0  3  1 ;
 2  4  0  6  2 ;
 0  0  0  1  5 ;
 0  0  0  2  10;
 0  0  0  3  15];

[m, cols] = size(An);
n = cols - 1; % So an so (n)
fprintf('Ma tran ban dau An (%d phuong trinh, %d an):\n', m, n);
disp(An);

i = 1; j = 1; buoc = 1;

% B2: Qua trinh thuan (Khu Gauss voi Partial Pivoting)
% Giu nguyen logic tim Pivot lon nhat de dam bao tinh chinh xac
while (i <= m) && (j <= n)
    [val, relative_p] = max(abs(An(i:m, j)));
    p = relative_p + i - 1;

    if val < 1e-10
        j = j + 1;
        continue;
    end

    if p ~= i
        An([i, p], :) = An([p, i], :);
    end

    for l = i + 1 : m
        heso = An(l, j) / An(i, j);
        An(l, :) = An(l, :) - heso * An(i, :);
    end

    % In tung buoc de theo doi
    fprintf('\n[BUOC %d] Khu tai vi tri (%d, %d)\n', buoc, i, j);
    disp(An);

    i = i + 1; j = j + 1;
    buoc = buoc + 1;
end

fprintf('============================================\n');
disp('MA TRAN SAU KHI KHU (DANG BAC THANG):');
disp(An);

% PHAN BIEN LUAN (B3 -> B7)

% B3 & B4: Khoi tao id[] va tim chi so cot khac 0 dau tien cua moi hang
id = zeros(1, m);
for r = 1:m
    idx = find(abs(An(r, :)) > 1e-10, 1);
    if ~isempty(idx)
        id(r) = idx;
    end
end

% B5: Tim imax (hang cuoi cung khac 0)
i_max = find(id ~= 0, 1, 'last');
if isempty(i_max); i_max = 0; end

% B6: Truong hop dac biet
if i_max == 0
    fprintf('B6: He co nghiem tuy y tren R^%d\n', n);
    return;
end

% B7: Bien luan
fprintf('B7: Bien luan nghiem (imax = %d, id(imax) = %d):\n', i_max, id(i_max));

if id(i_max) == n + 1
    % 7.1: Neu phan tu khac 0 dau tien nam o cot hang so
    fprintf('==> KET LUAN: HE VO NGHIEM (Xuat hien dang 0 = b).\n');
else
    % 7.2: He co nghiem
    k = n - i_max; % So an tu do

    % 7.2.2: Xac dinh cac an chinh va an tu do (pa)
    all_vars = 1:n;
    main_vars = id(1:i_max);
    pa = setdiff(all_vars, main_vars);

    if k == 0
        fprintf('==> KET LUAN: HE CO NGHIEM DUY NHAT.\n');
        % Qua trinh nghich (The nguoc)
        X = zeros(n, 1);
        for row = i_max:-1:1
            col_pivot = id(row);
            X(col_pivot) = (An(row, end) - An(row, col_pivot+1:n) * X(col_pivot+1:n)) / An(row, col_pivot);
        end
        disp('Nghiem X ='); disp(X);
    else
        fprintf('==> KET LUAN: HE CO VO SO NGHIEM (Phu thuoc %d tham so).\n', k);
        fprintf('Cac an tu do: '); disp(pa);
    end
end
