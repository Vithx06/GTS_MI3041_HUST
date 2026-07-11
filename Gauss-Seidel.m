% ==========================================
% CHUONG TRINH GIAI HE PHUONG TRINH GAUSS-SEIDEL
% ==========================================

clc; clear;

% ==========================================
% CAU HINH BAI TOAN (DAT TAI DAY)
% ==========================================
% Ma tran he so A
A = [
10 -1  2  0;
-1 11 -1  3;
2 -1 10 -1;
0  3 -1  8
];
% Ma tran/Vector B
B = [
6;
25;
-11;
15
];

% Xap xi ban dau X0
X0 = [0;
      0;
      0;
      0];

% Sai so muc tieu
epsilon = 1e-6;

% So buoc lap toi da
N0 = 100;

% Chon loai cheo troi (1: Troi hang, 2: Troi cot)
% Khong kiem tra lai dieu kien troi theo yeu cau
kieu_troi = 1;

% ==========================================
% BAT DAU THUAT TOAN
% ==========================================

n = size(A, 1);
m = size(B, 2);
I = eye(n);
T = diag(diag(A));
T_inv = inv(T);

fprintf('--- KHOI TAO THAM SO ---\n');

if kieu_troi == 1
    % Truong hop cheo troi hang
    C = I - T_inv * A;
    D = T_inv * B;

    mu = 0;
    for i = 1:n
        p_i = sum(abs(C(i, 1:i-1)));
        q_i = sum(abs(C(i, i:n)));
        mu_i = p_i / (1 - q_i);
        if mu_i > mu
            mu = mu_i;
        end
    end

    disp('Ma tran lap C =');
    disp(num2str(C, '%.7f    '));
    disp('Vector lap D =');
    disp(num2str(D, '%.7f    '));
    fprintf('He so co mu = %.7f\n\n', mu);

elseif kieu_troi == 2
    % Truong hop cheo troi cot
    M = I - T_inv * A;
    V = T_inv * B;
    C_col = (T - A) * T_inv;

    diag_A = abs(diag(A));
    lambda = max(diag_A) / min(diag_A);

    rho = 0;
    s_arr = zeros(n, 1);
    for j = 1:n
        t_j = sum(abs(M(1:j, j)));
        s_j = sum(abs(M(j+1:n, j)));
        s_arr(j) = s_j;
        rho_j = t_j / (1 - s_j);
        if rho_j > rho
            rho = rho_j;
        end
    end
    s = max(s_arr);

    disp('Ma tran lap M =');
    disp(num2str(M, '%.7f    '));
    disp('Vector lap V =');
    disp(num2str(V, '%.7f    '));
    fprintf('Tham so lambda = %.7f\n', lambda);
    fprintf('He so co rho = %.7f\n', rho);
    fprintf('He so s = %.7f\n\n', s);
end

% ==========================================
% VONG LAP GAUSS-SEIDEL
% ==========================================
fprintf('--- BAT DAU LAP ---\n');
X_k = X0;
k_final = 0;
eps_final = 0;

for k = 1:N0
    X_k_minus_1 = X_k;

    % Tinh toan X moi
    if kieu_troi == 1
        for i = 1:n
            sum1 = C(i, 1:i-1) * X_k(1:i-1, :);
            sum2 = C(i, i:n) * X_k_minus_1(i:n, :);
            X_k(i, :) = sum1 + sum2 + D(i, :);
        end
        err = (mu / (1 - mu)) * norm(X_k - X_k_minus_1, inf);
    else
        for i = 1:n
            sum1 = M(i, 1:i-1) * X_k(1:i-1, :);
            sum2 = M(i, i:n) * X_k_minus_1(i:n, :);
            X_k(i, :) = sum1 + sum2 + V(i, :);
        end
        err = (lambda * rho / ((1 - s) * (1 - rho))) * norm(X_k - X_k_minus_1, 1);
    end

    % In ket qua tung buoc
    fprintf('Buoc lap k = %d:\n', k);
    disp('Ma tran/Vector X =');
    disp(num2str(X_k, '%.7f    '));
    fprintf('Sai so eps = %.7f\n\n', err);

    k_final = k;
    eps_final = err;

    % Kiem tra dieu kien dung
    if err <= epsilon
        break;
    end
end

if k_final == N0 && eps_final > epsilon
    fprintf('Thong bao: Thuat toan khong hoi tu do vuot qua %d buoc lap.\n\n', N0);
end

% ==========================================
% TONG HOP KET QUA CUOI CUNG
% ==========================================
fprintf('==========================================\n');
fprintf('KET QUA CUOI CUNG\n');
fprintf('==========================================\n');
fprintf('So buoc lap da thuc hien: %d\n', k_final);
fprintf('Sai so dat duoc: %.7f\n', eps_final);
disp('Nghiem gan dung X tim duoc la:');
disp(num2str(X_k, '%.7f    '));
