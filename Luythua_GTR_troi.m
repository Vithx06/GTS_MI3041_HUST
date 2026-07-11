clc;
clear;

% ==================== CAU HINH BAI TOAN ====================
A = [4, 1, 1;
     1, 3, -1;
     1, -1, 2];          % Ma tran vuong A
Y0 = [1; 1; 1];          % Vector lap ban dau Y^(0)
epsilon = 1e-7;          % Sai so tieu chuan
N0 = 100;                % So buoc lap toi da
% ===========================================================

% Khoi tao bo nho cho lich su cac vector B
B_hist = zeros(size(A, 1), N0 + 1);
B_hist(:, 1) = Y0;
case_type = 3;           % Mac dinh la cap tri rieng phuc neu het vong lap ma khong thoat som
m_final = N0;

fprintf('BAT DAU CHUONG TRINH PHUONG PHAP LUY THUA\n\n');

for m = 1:N0
    % Tinh vector Y^(m+1)
    Y_next = A * B_hist(:, m);

    % Chuon hoa chuon vo cung de tim B^(m+1)
    max_val = max(abs(Y_next));
    B_next = Y_next / max_val;
    B_hist(:, m+1) = B_next;

    % In gia tri tim duoc o buoc lap hien tai
    fprintf('--- Buoc lap m = %d ---\n', m);
    fprintf('Vector Y^(%d):\n', m+1);
    for i = 1:length(Y_next)
        fprintf('  %.7f\n', Y_next(i));
    end
    fprintf('Vector B^(%d):\n', m+1);
    for i = 1:length(B_next)
        fprintf('  %.7f\n', B_next(i));
    end
    fprintf('\n');

    % Kiem tra hoi tu (voi m >= 3)
    if m >= 3
        err1 = max(abs(B_hist(:, m) - B_hist(:, m-1)));
        err2 = max(abs(B_hist(:, m) - B_hist(:, m-2)));

        if err1 <= epsilon
            case_type = 1;
            m_final = m;
            break;
        elseif err2 <= epsilon
            case_type = 2;
            m_final = m;
            break;
        end
    end
end

% ==================== TONG HOP KET QUA CUOI CUNG ====================
fprintf('==================================================\n');
fprintf('TONG HOP KET QUA CUOI CUNG\n');
fprintf('==================================================\n');
fprintf('So buoc lap da thuc hien: %d\n', m_final);

if case_type == 1
    fprintf('Truong hop: Tri rieng thuc troi\n');
    V = B_hist(:, m_final);
    [~, k] = max(abs(V));
    AV = A * V;
    lambda = AV(k) / V(k);
    V_norm = V / norm(V, 2);

    fprintf('Tri rieng lambda: %.7f\n', lambda);
    fprintf('Vector rieng da chuan hoa chuan 2:\n');
    for i = 1:length(V_norm)
        fprintf('  %.7f\n', V_norm(i));
    end

elseif case_type == 2
    fprintf('Truong hop: Hai tri rieng doi nhau\n');
    M = B_hist(:, m_final);
    AM = A * M;
    A2M = A * AM;
    [~, k] = max(abs(AM));

    lambda1 = sqrt(A2M(k) / M(k));
    lambda2 = -lambda1;

    V1 = AM + lambda1 * M;
    V2 = AM + lambda2 * M;

    V1_norm = V1 / norm(V1, 2);
    V2_norm = V2 / norm(V2, 2);

    fprintf('Tri rieng lambda1: %.7f\n', lambda1);
    fprintf('Vector rieng V1 tuong ung:\n');
    for i = 1:length(V1_norm)
        fprintf('  %.7f\n', V1_norm(i));
    end

    fprintf('Tri rieng lambda2: %.7f\n', lambda2);
    fprintf('Vector rieng V2 tuong ung:\n');
    for i = 1:length(V2_norm)
        fprintf('  %.7f\n', V2_norm(i));
    end

elseif case_type == 3
    fprintf('Truong hop: Cap tri rieng phuc lien hop\n');
    M = B_hist(:, m_final);
    AM = A * M;
    A2M = A * AM;

    % Lap he phuong trinh tuyen tinh tu dong 1 va dong 2
    Sys_Mat = [AM(1), M(1); AM(2), M(2)];
    Sys_Vec = [-A2M(1); -A2M(2)];
    bc = Sys_Mat \ Sys_Vec;
    b = bc(1); c = bc(2);

    % Giai phuong trinh bac hai de tim tri rieng phuc
    roots_lambda = roots([1, b, c]);
    lambda1 = roots_lambda(1);
    lambda2 = roots_lambda(2);

    V1 = AM - lambda1 * M;
    V2 = AM - lambda2 * M;

    V1_norm = V1 / norm(V1, 2);
    V2_norm = V2 / norm(V2, 2);

    fprintf('Tri rieng lambda1: %.7f %+.7fi\n', real(lambda1), imag(lambda1));
    fprintf('Vector rieng V1 tuong ung:\n');
    for i = 1:length(V1_norm)
        fprintf('  %.7f %+.7fi\n', real(V1_norm(i)), imag(V1_norm(i)));
    end

    fprintf('Tri rieng lambda2: %.7f %+.7fi\n', real(lambda2), imag(lambda2));
    fprintf('Vector rieng V2 tuong ung:\n');
    for i = 1:length(V2_norm)
        fprintf('  %.7f %+.7fi\n', real(V2_norm(i)), imag(V2_norm(i)));
    end
end
