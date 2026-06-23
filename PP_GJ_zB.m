% =========================================================================
% CHUONG TRINH GIAI PHUONG TRINH MA TRAN AX = B BANG GAUSS-JORDAN
% (Phan tu troi toan phan - Khong hoan doi hang - B co z cot)
% =========================================================================
clc; clear;

% B1: Cau hinh du lieu
A = [-3,5,5,-11,0;
      -2,11,17,-35,2;
      -4,-9,2,6,12 ;
      9,-3,25,-39,16;
     6,4,1,-7,-6 ];

B = [ -4   -1 ;
      -7   6 ;
     7   -8 ;
     8  19 ;
     -2 10];

A_bar = [A, B]; % Ma tran bo sung [A | B]

[m, cols] = size(A_bar);
z = size(B, 2); % So cot cua ma tran B
n = cols - z;   % So an (So cot cua A)

fprintf('Ma tran bo sung ban dau (%d phuong trinh, ma tran A co %d cot, B co %d cot):\n', m, n, z);
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
    % Xac dinh cac hang va cot cua A chua bi khu
    valid_rows = setdiff(1:m, usedR);
    valid_cols = setdiff(1:n, usedC); % Chi tim Pivot trong vung ma tran A

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
    % Phep chia tu dong ap dung cho ca z cot cua B
    A_bar(p, :) = A_bar(p, :) / A_bar(p, q);

    % 2.3.3: Khu cot q (Tren va duoi Pivot)
    for l = 1:m
        if l ~= p && abs(A_bar(l, q)) > 1e-10
            heso = A_bar(l, q);
            % Phep tru tu dong ap dung cho toan bo hang, bao gom ca B
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
    % Kiem tra xem co bat ky cot nao trong B xuat hien 0 = b != 0 khong
    if any(abs(A_bar(r, n+1:end)) > 1e-10)
        vo_nghiem = true;
        break;
    end
end

if vo_nghiem
    fprintf('==> KET LUAN: PHUONG TRINH MA TRAN VO NGHIEM.\n');
    fprintf('(Xuat hien mâu thuan 0 = b tai it nhat 1 cot cua B).\n');
else
    % B4: Xac dinh hang
    s = length(usedR);
    if s == 0
        fprintf('==> KET LUAN: He co ma tran nghiem X tuy y tren khong gian R^(%dx%d)\n', n, z);
        return;
    end

    % B5: Phan loai an so
    k = n - s;
    pa = setdiff(1:n, usedC); % Cac an tu do

    % B6: Ket luan nghiem
    if k == 0
        fprintf('==> KET LUAN: CO MA TRAN NGHIEM DUY NHAT.\n');
        X = zeros(n, z); % Khoi tao ma tran nghiem
        for p = usedR
            % Doc truc tiep tat ca z cot tu A_bar
            X(id(p), :) = A_bar(p, n+1:end);
        end

        disp('Ma tran nghiem X = ');
        disp(X);
    else
        fprintf('==> KET LUAN: CO VO SO MA TRAN NGHIEM (Moi cot phu thuoc %d tham so).\n', k);
        fprintf('Danh sach cac hang cua X dong vai tro an tu do: '); disp(pa);

        % Sap xep lai the hien theo thu tu an chinh tu nho den lon cho dep mat
        sorted_id = sort(id(usedR));

        % In bieu dien nghiem cho tung cot cua ma tran X
        for c = 1:z
            fprintf('\n--- Bieu dien nghiem cho COT %d cua ma tran X ---\n', c);
            for q = sorted_id
                p = find(id == q); % Tim lai hang tuong ung voi an chinh q

                % In he so tu do cua cot c hien tai
                fprintf('  X(%d,%d) = %10.6f ', q, c, A_bar(p, n+c));

                % Ghep cac an tu do
                for t = pa
                    heso_tu_do = -A_bar(p, t);
                    if abs(heso_tu_do) > 1e-10
                        if heso_tu_do > 0
                            fprintf('+ %.6f*X(%d,%d) ', heso_tu_do, t, c);
                        else
                            fprintf('- %.6f*X(%d,%d) ', abs(heso_tu_do), t, c);
                        end
                    end
                end
                fprintf('\n');
            end
        end
    end
end
