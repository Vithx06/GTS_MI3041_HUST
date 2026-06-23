% =========================================================================
% CHUONG TRINH GIAI HE PHUONG TRINH BANG GAUSS-JORDAN
% (Phan tu troi toan phan - Khong hoan doi hang)
% =========================================================================
clc; clear;

% B1: Cau hinh du lieu (Ma tran bo sung A_bar = [A | B])
% Vi du: He 4 phuong trinh, 5 an
A_bar = [ 0   2  -1   3  2 ;
          2   5  -5   8  2 ;
          5   4   4   6  13 ;
          4   8  -2  14  14 ];

[m, cols] = size(A_bar);
n = cols - 1; % So an so

fprintf('Ma tran bo sung ban dau (%d phuong trinh, %d an):\n', m, n);
disp(A_bar);

% Khoi tao bien (Buoc 1)
usedR = []; % Tap cac hang da khu
usedC = []; % Tap cac cot da khu
id = zeros(1, m); % Mang luu vet an chinh cua tung hang
buoc = 1;

% =========================================================================
% B2: QUA TRINH KHU GAUSS-JORDAN
% =========================================================================
while length(usedR) < m && length(usedC) < n
    % Xac dinh cac hang va cot chua bi khu
    valid_rows = setdiff(1:m, usedR);
    valid_cols = setdiff(1:n, usedC);

    % Trich xuat vung ma tran chua khu de tim Pivot
    subA = A_bar(valid_rows, valid_cols);

    % 2.1: Tim phan tu khu (Pivot)
    val_max = max(abs(subA(:)));

    % 2.2: Kiem tra dung
    if val_max < 1e-10
        break; % Tat ca he so con lai deu bang 0 -> Xong qua trinh khu
    end

    % Uu tien 1: Tim phan tu co tri tuyet doi bang 1 (Sai so 1e-10)
    [r_ones, c_ones] = find(abs(abs(subA) - 1) < 1e-10);

    if ~isempty(r_ones)
        p_rel = r_ones(1);
        q_rel = c_ones(1);
        loai_uu_tien = 1;
    else
        % Uu tien 2: Chon phan tu lon nhat
        [~, loc_max] = max(abs(subA(:)));
        [p_rel, q_rel] = ind2sub(size(subA), loc_max);
        loai_uu_tien = 2;
    end

    % Chuyen doi toa do tuong doi (subA) thanh toa do tuyet doi (A_bar)
    p = valid_rows(p_rel);
    q = valid_cols(q_rel);

    fprintf('\n[BUOC %d] -------------------------------------------\n', buoc);
    fprintf('Tim thay Pivot = %.4f tai (hang %d, cot %d) - Uu tien %d\n', A_bar(p, q), p, q, loai_uu_tien);

    % 2.3.1: Ghi vet
    usedR = [usedR, p];
    usedC = [usedC, q];
    id(p) = q;

    % 2.3.2: Chuan hoa hang p
    A_bar(p, :) = A_bar(p, :) / A_bar(p, q);

    % 2.3.3: Khu cot q (Tren va duoi Pivot)
    for l = 1:m
        if l ~= p && abs(A_bar(l, q)) > 1e-10
            heso = A_bar(l, q);
            A_bar(l, :) = A_bar(l, :) - heso * A_bar(p, :);
        end
    end

    disp('Ma tran sau khi chuan hoa va khu:');
    disp(A_bar);
    buoc = buoc + 1;
end

% =========================================================================
% BIEN LUAN KET QUA (B3 -> B6)
% =========================================================================
fprintf('============================================\n');
disp('MA TRAN KET QUA GAUSS-JORDAN:');
disp(A_bar);

% B3: Kiem tra vo nghiem
vo_nghiem = false;
unusedR = setdiff(1:m, usedR); % Cac hang chua bi chon lam Pivot
for r = unusedR
    if abs(A_bar(r, n+1)) > 1e-10
        vo_nghiem = true;
        break;
    end
end

if vo_nghiem
    fprintf('==> KET LUAN: HE VO NGHIEM (Xuat hien phuong trinh 0 = b != 0).\n');
else
    % B4: Xac dinh hang
    s = length(usedR);
    if s == 0
        fprintf('==> KET LUAN: He co nghiem tuy y tren R^%d\n', n);
        return;
    end

    % B5: Phan loai an so
    k = n - s;
    pa = setdiff(1:n, usedC); % Cac an tu do

    % B6: Ket luan nghiem
    if k == 0
        fprintf('==> KET LUAN: HE CO NGHIEM DUY NHAT.\n');
        X = zeros(n, 1);
        for p = usedR
            X(id(p)) = A_bar(p, n+1);
        end

        fprintf('Nghiem X = \n');
        for v = 1:n
            fprintf('  x_%d = %10.6f\n', v, X(v));
        end
    else
        fprintf('==> KET LUAN: HE CO VO SO NGHIEM (Phu thuoc %d tham so).\n', k);
        fprintf('Danh sach cac an tu do: '); disp(pa);

        fprintf('\nBieu dien cac an chinh theo an tu do:\n');
        % Sap xep lai the hien theo thu tu an chinh tu nho den lon cho dep mat
        sorted_id = sort(id(usedR));
        for q = sorted_id
            p = find(id == q); % Tim lai hang tuong ung voi an chinh q
            fprintf('  x_%d = %10.6f ', q, A_bar(p, n+1));
            for t = pa
                heso_tu_do = -A_bar(p, t);
                if abs(heso_tu_do) > 1e-10
                    if heso_tu_do > 0
                        fprintf('+ %.6f*x_%d ', heso_tu_do, t);
                    else
                        fprintf('- %.6f*x_%d ', abs(heso_tu_do), t);
                    end
                end
            end
            fprintf('\n');
        end
    end
end
