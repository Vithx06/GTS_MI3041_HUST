% Xoa man hinh va du lieu cu
clc;
clear all;

% =========================================================================
% CAU HINH BAI TOAN (INPUT)
% Thay doi gia tri ma tran A va B tai day
% =========================================================================
A = [1, 2, 3;
     2, 5, 2;
     3, 1, 5];

B = [14;
     18;
     20];

% Ham ho tro in ma tran voi 7 chu so sau dau phay
function in_matran(ten_matran, M)
    fprintf('%s =\n', ten_matran);
    [dong, cot] = size(M);
    for i = 1:dong
        for j = 1:cot
            fprintf('%14.7f ', M(i,j));
        end
        fprintf('\n');
    end
    fprintf('\n');
end

% =========================================================================
% THUC THI THUAT TOAN
% =========================================================================

% --- BUOC 1: KIEM TRA KICH THUOC ---
fprintf('--- BUOC 1: KIEM TRA KICH THUOC ---\n');
[m, n] = size(A);
[mb, pb] = size(B);
if m ~= mb
    fprintf('Loi: Kich thuoc ma tran khong tuong thich (so hang cua A va B phai bang nhau).\n');
    return;
else
    fprintf('Kich thuoc hop le: m = %d, n = %d, p = %d\n\n', m, n, pb);
end

% --- BUOC 2: DUA VE HE VUONG ---
fprintf('--- BUOC 2: DUA VE HE VUONG ---\n');
if m > n
    fprintf('Trang thai: m > n. Dua ve he tuong duong: A^T A X = A^T B\n');
    A_vuong = A' * A;
    B_vuong = A' * B;
elseif m < n
    fprintf('Trang thai: m < n. Dua ve he tuong duong: (A A^T) Y = B\n');
    A_vuong = A * A';
    B_vuong = B;
else
    fprintf('Trang thai: m = n. Giu nguyen he phuong trinh.\n');
    A_vuong = A;
    B_vuong = B;
end
N = size(A_vuong, 1);
in_matran('Ma tran he vuong dang xet', A_vuong);

% --- BUOC 3: KHOI TAO PHAN RA PA = LU ---
fprintf('--- BUOC 3: KHOI TAO ---\n');
L = eye(N);
U = eye(N);
P = eye(N);
fprintf('Hoan tat khoi tao L, U, P la ma tran don vi cap %d.\n\n', N);

% --- BUOC 4: KIEM TRA DINH THUC CON CHINH ---
% Da vo hieu hoa de tranh loi sai so thap phan.
% Thuat toan se dua vao kiem tra phan tu troi o Buoc 5.2.
% fprintf('--- BUOC 4: KIEM TRA DINH THUC CON CHINH (Da bo qua) ---\n\n');

% --- BUOC 5: VONG LAP CHINH PHAN RA ---
fprintf('--- BUOC 5: VONG LAP CHINH PHAN RA ---\n');
for k = 1:N
    fprintf('>> Vong lap k = %d:\n', k);

    % 5.1 Tinh cot k cua L
    for i = k:N
        tong_LU = 0;
        for t = 1:k-1
            tong_LU = tong_LU + L(i,t) * U(t,k);
        end
        L(i,k) = A_vuong(i,k) - tong_LU;
        fprintf('   Tinh L(%d,%d) = %.7f\n', i, k, L(i,k));
    end

    % 5.2 Chon phan tu troi
    [gia_tri_max, vi_tri_max] = max(abs(L(k:N, k)));
    p = vi_tri_max + k - 1;
    fprintf('   Phan tu troi nam o hang p = %d (gia tri: %.7f)\n', p, L(p,k));

    if L(p,k) == 0
        fprintf('Loi: Ma tran suy bien, thuat toan ket thuc.\n');
        return;
    end

    % 5.3 Hoan vi hang
    if p ~= k
        fprintf('   Thuc hien hoan vi hang %d va %d.\n', k, p);

        % Doi cho A_vuong
        tam = A_vuong(k,:); A_vuong(k,:) = A_vuong(p,:); A_vuong(p,:) = tam;

        % Doi cho P
        tam = P(k,:); P(k,:) = P(p,:); P(p,:) = tam;

        % Doi cho L (cac phan tu tu cot 1 den k)
        tam = L(k, 1:k); L(k, 1:k) = L(p, 1:k); L(p, 1:k) = tam;
    end

    % 5.4 Tinh hang k cua U
    for j = k:N
        tong_LU = 0;
        for t = 1:k-1
            tong_LU = tong_LU + L(k,t) * U(t,j);
        end
        U(k,j) = (A_vuong(k,j) - tong_LU) / L(k,k);
        if j > k
            fprintf('   Tinh U(%d,%d) = %.7f\n', k, j, U(k,j));
        end
    end

    % 5.5 In trang thai L va U hien tai
    fprintf('\n   --- Trang thai ma tran L va U sau vong lap k = %d ---\n', k);
    in_matran('   L', L);
    in_matran('   U', U);
    fprintf('------------------------------------------------------\n\n');
end

% --- BUOC 6: GIAI HE PHUONG TRINH ---
fprintf('--- BUOC 6: GIAI HE ---\n');
B_phay = P * B_vuong;
in_matran('B_phay (P * B_vuong)', B_phay);

Y = zeros(N, pb);
fprintf('The tien tinh ma tran Y:\n');
for col = 1:pb
    for i = 1:N
        tong_Ly = 0;
        for j = 1:i-1
            tong_Ly = tong_Ly + L(i,j) * Y(j, col);
        end
        Y(i, col) = (B_phay(i, col) - tong_Ly) / L(i,i);
        fprintf('   Y(%d,%d) = %.7f\n', i, col, Y(i, col));
    end
end
fprintf('\n');

Z = zeros(N, pb); % Day la ma tran nghiem cua he vuong (tuong duong X hoac Y tuy truong hop o Buoc 2)
fprintf('The nguoc tinh ma tran nghiem he vuong:\n');
for col = 1:pb
    for i = N:-1:1
        tong_Uz = 0;
        for j = i+1:N
            tong_Uz = tong_Uz + U(i,j) * Z(j, col);
        end
        Z(i, col) = Y(i, col) - tong_Uz;
        fprintf('   Nghiem_vuong(%d,%d) = %.7f\n', i, col, Z(i, col));
    end
end
fprintf('\n');

% --- BUOC 7: CHUYEN DOI NGHIEM ---
fprintf('--- BUOC 7: CHUYEN DOI NGHIEM ---\n');
if m < n
    X_kq = A' * Z;
    fprintf('Truong hop m < n, tinh nghiem thuc te X = A^T * Y_he_vuong\n');
else
    X_kq = Z;
    fprintf('Khong thuoc truong hop m < n. Nghiem giu nguyen.\n');
end
fprintf('\n');

% --- BUOC 8: TRA VE KET QUA TONG HOP ---
fprintf('=========================================================================\n');
fprintf('--- BUOC 8: TONG HOP CAC SO LIEU QUAN TRONG ---\n');
fprintf('=========================================================================\n');
in_matran('MA TRAN HOAN VI (P)', P);
in_matran('MA TRAN TAM GIAC DUOI (L)', L);
in_matran('MA TRAN TAM GIAC TREN (U)', U);
in_matran('MA TRAN NGHIEM (X)', X_kq);
