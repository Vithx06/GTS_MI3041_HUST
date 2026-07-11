clc; clear all;

% ====================================================
% CAU HINH BAI TOAN
% ====================================================
A =   [
    33 -1 2 1 -3 2;
    -5 33 -1 5 2 -3;
    1 3 33 2 2 1;
    2 1 -1 33 0 4;
    3 -8 4 3 33 2;
    2 3 1 4 2 33
];
B =  [
    29.606 24.8297;
    3.5764 8.5379;
    26.9654 11.342;
    34.3139 31.1558;
    16.2633 15.5075;
    19.7982 29.9398
];
X0 = [0 0;
      0 0;
      0 0;
      0 0;
      0 0;
      0 0];
epsilon = 1e-6;
N0 = 100;

% Chon loai cheo troi muon su dung:
% 1: Cheo troi hang
% 2: Cheo troi cot
lua_chon_troi = 1;

% ====================================================
% THUAT TOAN JACOBI
% ====================================================

% --- B1: KHOI TAO ---
fprintf('--- B1: KHOI TAO ---\n');
fprintf('Ma tran A:\n'); disp(num2str(A, '%12.7f'));
fprintf('Vecto B:\n'); disp(num2str(B, '%12.7f'));
fprintf('Xap xi ban dau X0:\n'); disp(num2str(X0, '%12.7f'));
fprintf('Sai so epsilon = %.7f\n', epsilon);
fprintf('So buoc lap N0 = %d\n', N0);

if lua_chon_troi == 1
    fprintf('Che do dang chon: CHEO TROI HANG\n\n');
elseif lua_chon_troi == 2
    fprintf('Che do dang chon: CHEO TROI COT\n\n');
else
    fprintf('Loi: lua_chon_troi chi nhan gia tri 1 hoac 2. Ket thuc.\n');
    return;
end

n = size(A, 1);
T = diag(diag(A));
fprintf('Ma tran T:\n'); disp(num2str(T, '%12.7f'));
fprintf('\n');

% --- B2: THIET LAP CAC MA TRAN LUU CHUYEN ---
fprintf('--- B2: THIET LAP MA TRAN ---\n');
invT = inv(T);

if lua_chon_troi == 1
    C = eye(n) - invT * A;
    D = invT * B;
    q = norm(C, inf);

    fprintf('Ma tran C:\n'); disp(num2str(C, '%12.7f'));
    fprintf('Vecto D:\n'); disp(num2str(D, '%12.7f'));
    fprintf('He so co q = %.7f\n', q);

    if q >= 1 || q == 0
        fprintf('Khong thoa man dieu kien he so co q. Ket thuc.\n');
        return;
    end

elseif lua_chon_troi == 2
    M = eye(n) - invT * A;
    V = invT * B;
    C_cot = (T - A) * invT;

    duong_cheo = abs(diag(A));
    lambda = max(duong_cheo) / min(duong_cheo);
    q = norm(C_cot, 1);

    fprintf('Ma tran M:\n'); disp(num2str(M, '%12.7f'));
    fprintf('Vecto V:\n'); disp(num2str(V, '%12.7f'));
    fprintf('Ma tran C (theo cot):\n'); disp(num2str(C_cot, '%12.7f'));
    fprintf('Tham so lambda = %.7f\n', lambda);
    fprintf('He so co q = %.7f\n', q);

    if q >= 1 || q == 0
        fprintf('Khong thoa man dieu kien he so co q. Ket thuc.\n');
        return;
    end
end
fprintf('\n');

% --- B3 & B4: QUA TRINH LAP ---
fprintf('--- B3 & B4: QUA TRINH LAP ---\n');
k = 1; % B3
X_truoc = X0;
hoi_tu = 0;

while k <= N0 % B4
    % 4.1: Tinh xap xi
    if lua_chon_troi == 1
        X_hien_tai = C * X_truoc + D;
        % 4.2: Tinh sai so
        eps = (q / (1 - q)) * norm(X_hien_tai - X_truoc, inf);
    else
        X_hien_tai = M * X_truoc + V;
        % 4.2: Tinh sai so
        eps = (lambda * q / (1 - q)) * norm(X_hien_tai - X_truoc, 1);
    end

    % In ket qua buoc lap
    fprintf('Buoc k = %d:\n', k);
    fprintf('X:\n'); disp(num2str(X_hien_tai, '%12.7f'));
    fprintf('Sai so eps = %.7f\n\n', eps);

    % 4.3: Kiem tra dung
    if eps <= epsilon
        hoi_tu = 1;
        break; % Chuyen sang B6
    end

    % 4.4: Cap nhat
    X_truoc = X_hien_tai;
    k = k + 1;
end

% --- B5 & B6: TONG HOP KET QUA ---
fprintf('--- TONG HOP KET QUA ---\n');
if hoi_tu == 1
    % B6
    fprintf('Trang thai: HOI TU\n');
    fprintf('So buoc lap thuc hien k = %d\n', k);
    fprintf('Nghiem gan dung X:\n'); disp(num2str(X_hien_tai, '%12.7f'));
    fprintf('Sai so dat duoc eps = %.7f\n', eps);
else
    % B5
    fprintf('Trang thai: KHONG HOI TU do vuot qua so buoc lap N0 = %d.\n', N0);
end
