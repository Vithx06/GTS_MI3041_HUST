% ==========================================
% CAU HINH BAI TOAN (Dau vao)
% ==========================================
clc; clear;

A = A = [2, 1;
     1, -3];
B = [5;
     -1];

function in_ma_tran(M, ten_ma_tran)
    fprintf('%s :\n', ten_ma_tran);
    [r, c] = size(M);
    for i = 1:r
        for j = 1:c
            fprintf('%12.7f ', M(i,j));
        end
        fprintf('\n');
    end
    fprintf('\n');
end

fprintf('CHUONG TRINH GIAI HE PHUONG TRINH BANG PHAN RA CHOLESKY\n');
fprintf('========================================================\n\n');

% ==========================================
% B1: KHOI TAO
% ==========================================
fprintf('--- B1: KHOI TAO ---\n');
[m, n] = size(A);
[mB, p] = size(B);

if m ~= mB
    fprintf('Loi: So hang cua A va B phai bang nhau.\n');
    return;
end

in_ma_tran(A, 'Ma tran A');
in_ma_tran(B, 'Ma tran B');

% ==========================================
% B2: DUA VE HE VUONG
% ==========================================
fprintf('--- B2: DUA VE HE VUONG ---\n');
if m > n
    fprintf('Truong hop m > n. Giai he (A^T * A)X = A^T * B\n');
    A_prime = A' * A;
    B_prime = A' * B;
elseif m < n
    fprintf('Truong hop m < n. Giai he (A * A^T)Y = B\n');
    A_prime = A * A';
    B_prime = B;
else
    fprintf('Truong hop m = n.\n');
    A_prime = A;
    B_prime = B;
end

in_ma_tran(A_prime, 'Ma tran A phay');
in_ma_tran(B_prime, 'Ma tran B phay');

% ==========================================
% B3: KIEM TRA KHA NGHICH
% ==========================================
fprintf('--- B3: KIEM TRA KHA NGHICH ---\n');
det_A_prime = det(A_prime);
fprintf('Dinh thuc cua A phay = %.7f\n', det_A_prime);

if abs(det_A_prime) < 1e-12
    fprintf('Ma tran khong kha nghich, dung chuong trinh.\n');
    return;
end

if m >= n
    N = n;
else
    N = m;
end
fprintf('Cap cua ma tran vuong A phay, N = %d\n\n', N);

% ==========================================
% B4: PHAN RA CHOLESKY MA TRAN A'
% ==========================================
fprintf('--- B4: PHAN RA CHOLESKY ---\n');
Q = zeros(N, N);

for i = 1:N
    fprintf('>> Vong lap hang i = %d:\n', i);

    % 4.1 Tinh phan tu duong cheo
    s_diag = 0;
    for k = 1:i-1
        s_diag = s_diag + Q(k,i)^2;
    end

    if A_prime(i,i) - s_diag <= 0
        fprintf('Ma tran khong xac dinh duong. Dung thuat toan.\n');
        return;
    end

    Q(i,i) = sqrt(A_prime(i,i) - s_diag);
    fprintf('   Phan tu duong cheo Q(%d,%d) = %.7f\n', i, i, Q(i,i));

    % 4.2 Tinh cac phan tu con lai tren hang i
    for j = i+1:N
        s_row = 0;
        for k = 1:i-1
            s_row = s_row + Q(k,i) * Q(k,j);
        end
        Q(i,j) = (A_prime(i,j) - s_row) / Q(i,i);
        fprintf('   Phan tu Q(%d,%d) = %.7f\n', i, j, Q(i,j));
    end
    fprintf('\n');
end

in_ma_tran(Q, 'Ma tran tam giac tren Q sau B4');

% ==========================================
% B5: GIAI HE TRUNG GIAN VA TIM NGHIEM
% ==========================================
fprintf('--- B5: GIAI HE TRUNG GIAN ---\n');
Z = zeros(N, p);

fprintf('5.1 Giai he tam giac duoi Q^T * Z = B phay (Phuong phap the tien):\n');
for c = 1:p
    fprintf(' >> Cot c = %d:\n', c);
    for i = 1:N
        sum_val = 0;
        for k = 1:i-1
            sum_val = sum_val + Q(k,i) * Z(k,c);
        end
        Z(i,c) = (B_prime(i,c) - sum_val) / Q(i,i);
        fprintf('    Z(%d,%d) = %.7f\n', i, c, Z(i,c));
    end
end
fprintf('\n');
in_ma_tran(Z, 'Ma tran trung gian Z hoan chinh');

fprintf('5.2 Giai he tam giac tren (Phuong phap the nguoc):\n');
if m >= n
    X = zeros(N, p);
    for c = 1:p
        fprintf(' >> Cot c = %d:\n', c);
        for i = N:-1:1
            sum_val = 0;
            for k = i+1:N
                sum_val = sum_val + Q(i,k) * X(k,c);
            end
            X(i,c) = (Z(i,c) - sum_val) / Q(i,i);
            fprintf('    X(%d,%d) = %.7f\n', i, c, X(i,c));
        end
    end
else
    Y = zeros(N, p);
    for c = 1:p
        fprintf(' >> Cot c = %d:\n', c);
        for i = N:-1:1
            sum_val = 0;
            for k = i+1:N
                sum_val = sum_val + Q(i,k) * Y(k,c);
            end
            Y(i,c) = (Z(i,c) - sum_val) / Q(i,i);
            fprintf('    Y(%d,%d) = %.7f\n', i, c, Y(i,c));
        end
    end
    fprintf('\n');
    in_ma_tran(Y, 'Ma tran trung gian Y hoan chinh');
    X = A' * Y;
end
fprintf('\n');

% ==========================================
% B6: TRA VE NGHIEM
% ==========================================
fprintf('--- B6: KET QUA NGHIEM ---\n');
in_ma_tran(X, 'Ma tran nghiem X');

% ==========================================
% TONG HOP SO LIEU QUAN TRONG
% ==========================================
fprintf('========================================================\n');
fprintf('TONG HOP SO LIEU QUAN TRONG:\n');
fprintf('- Kich thuoc ma tran A ban dau: %d hang x %d cot\n', m, n);
fprintf('- Dinh thuc he so thu gon det(A phay): %.7f\n', det_A_prime);
fprintf('- Thuat toan su dung: Phan ra Cholesky\n\n');
in_ma_tran(Q, '- Ma tran phan ra Q (Tam giac tren)');
in_ma_tran(X, '- Nghiem cuoi cung he AX = B (Ma tran X)');
fprintf('========================================================\n');
