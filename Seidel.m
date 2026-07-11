% =========================================================================
% CHUONG TRINH GIAI HE PHUONG TRINH BANG PHUONG PHAP LAP SEIDEL
% =========================================================================

clc; clear;

% --- CAU HINH BAI TOAN (DAT TAI DAY) ---
% Ma tran A (Kich thuoc n x n)
A = [
    0.9600  -0.0700   0.0100   0        0.0500;
    0.1000   0.9600   0.0200   0.0100  -0.0400;
    0.0500   0.0400   0.9400  -0.0300  -0.0300;
    0.1000  -0.0900  -0.0600   0.9600   0.0700;
    0.0800   0.1000   0.0700  -0.0500   1.0800
];

% Ma tran/Vector B (Kich thuoc n x m)
B =  [
   -7;
    6;
   -4;
    1;
   -7
];

% Xap xi dau X0 (Cung kich thuoc voi B)
X0 = [0; 0; 0; 0; 0];

% Sai so cho phep
epsilon = 1e-5;

% So buoc lap toi da
N0 = 5;

% Lua chon chuan: 1 (Cho chuan 1) hoac inf (Cho chuan vo cung)
chuan = inf;
% ---------------------------------------

% B1: KHOI TAO
[n, m] = size(B);
I = eye(n);
C = I - A;

disp('======================================================');
disp('      PHUONG PHAP LAP SEIDEL GIAI HE PHUONG TRINH     ');
disp('======================================================');
disp('--- BUOC 1: KHOI TAO ---');
disp('Ma tran C = I - A:');
for r = 1:n
    for c = 1:n
        fprintf('%12.7f ', C(r,c));
    end
    fprintf('\n');
end
fprintf('\n');

% B2: TINH TOAN HE SO VA KIEM TRA HOI TU
disp('--- BUOC 2: KIEM TRA DIEU KIEN TRANG THAI ---');
if chuan == inf
    normC = norm(C, inf);
    fprintf('Chuan vo cung cua C = %.7f\n', normC);

    if normC >= 1
        disp('Loi: Thuat toan co the khong hoi tu vi chuan C >= 1. KET THUC.');
        return;
    end

    p = zeros(n, 1);
    q = zeros(n, 1);
    for i = 1:n
        p(i) = sum(abs(C(i, 1:i-1)));
        q(i) = sum(abs(C(i, i:n)));
    end
    mu = max(p ./ (1 - q));
    fprintf('He so co mu = %.7f\n\n', mu);

elseif chuan == 1
    normC = norm(C, 1);
    fprintf('Chuan 1 cua C = %.7f\n', normC);

    if normC >= 1
        disp('Loi: Thuat toan co the khong hoi tu vi chuan C >= 1. KET THUC.');
        return;
    end

    t = zeros(n, 1);
    s_vec = zeros(n, 1);
    for j = 1:n
        t(j) = sum(abs(C(1:j, j)));
        if j < n
            s_vec(j) = sum(abs(C(j+1:n, j)));
        else
            s_vec(j) = 0;
        end
    end
    s = max(s_vec);
    rho = max(t ./ (1 - s_vec));
    fprintf('He so s = %.7f\n', s);
    fprintf('He so co rho = %.7f\n\n', rho);
else
    disp('Loi: Chuan khong hop le. Vui long chon 1 hoac inf tai muc cau hinh.');
    return;
end

% B3: KHOI TAO BUOC LAP
k = 1;
X_prev = X0;
X_curr = zeros(n, m);

disp('--- BUOC 3 & 4: QUA TRINH LAP ---');
% B4: VONG LAP
while k <= N0
    fprintf('\n--- Vong lap k = %d ---\n', k);

    % 4.1: Tinh cac hang cua X(k)
    for i = 1:n
        sum1 = zeros(1, m);
        sum2 = zeros(1, m);

        if i > 1
            sum1 = C(i, 1:i-1) * X_curr(1:i-1, :);
        end
        sum2 = C(i, i:n) * X_prev(i:n, :);

        X_curr(i, :) = sum1 + sum2 + B(i, :);
    end

    % In ma tran X dat duoc o buoc k
    fprintf('Ma tran/Vector X(%d):\n', k);
    for r = 1:n
        for c = 1:m
            fprintf('%12.7f ', X_curr(r,c));
        end
        fprintf('\n');
    end

    % 4.2: Tinh sai so eps
    if chuan == inf
        diff_norm = norm(X_curr - X_prev, inf);
        eps_val = (mu / (1 - mu)) * diff_norm;
    else
        diff_norm = norm(X_curr - X_prev, 1);
        eps_val = (rho / ((1 - s) * (1 - rho))) * diff_norm;
    end

    fprintf('Sai so eps = %.7f\n', eps_val);

    % 4.3: Kiem tra dieu kien dung
    if eps_val <= epsilon
        disp('-> Sai so da thoa man dieu kien. Dung lap.');
        break;
    end

    % 4.4: Tang k
    X_prev = X_curr;
    k = k + 1;
end

% B5 & B6: KET LUAN
disp('\n======================================================');
if k > N0
    disp('--- BUOC 5: THONG BAO ---');
    fprintf('Thuat toan KHONG HOI TU do vuot qua %d buoc lap.\n', N0);
else
    disp('--- BUOC 6: KET QUA TONG HOP ---');
    fprintf('So buoc lap thuc hien (k) : %d\n', k);
    fprintf('Sai so dat duoc (eps)     : %.7f\n', eps_val);
    disp('Nghiem gan dung X tim duoc:');
    for r = 1:n
        for c = 1:m
            fprintf('%12.7f ', X_curr(r,c));
        end
        fprintf('\n');
    end
end
disp('======================================================');
