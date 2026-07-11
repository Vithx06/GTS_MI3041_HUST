clc;
clear;

% =========================================================================
% 1. DINH NGHIA CAC HAM BO TRO
% =========================================================================
function disp_matrix(M)
    [r, c] = size(M);
    for row = 1:r
        for col = 1:c
            val = M(row, col);
            if imag(val) == 0
                fprintf('%14.7f ', real(val));
            elseif imag(val) >= 0
                fprintf('%14.7f +%11.7fi ', real(val), imag(val));
            else
                fprintf('%14.7f -%11.7fi ', real(val), abs(imag(val)));
            end
        end
        fprintf('\n');
    end
end

function disp_vector(V)
    for idx = 1:length(V)
        val = V(idx);
        if imag(val) == 0
            fprintf('%14.7f\n', real(val));
        elseif imag(val) >= 0
            fprintf('%14.7f +%11.7fi\n', real(val), imag(val));
        else
            fprintf('%14.7f -%11.7fi\n', real(val), abs(imag(val)));
        end
    end
end

function A_prime = Xuong_Thang(A_mat, V_vec)
    n_dim = size(A_mat, 1);
    [~, k_idx] = max(abs(V_vec));
    Theta = eye(n_dim);
    Theta(:, k_idx) = Theta(:, k_idx) - V_vec;
    Theta(k_idx, k_idx) = 0;
    A_prime = Theta * A_mat;
end

% =========================================================================
% 2. CAU HINH BAI TOAN
% =========================================================================
A_init = [4, -1, 1;
          -1, 3, -2;
          1, -2, 3];

Y0 = [1; 1; 1];
epsilon = 1e-7;
N0 = 50;

% =========================================================================
% 3. THUAT TOAN CHINH
% =========================================================================
A = A_init;
n = size(A, 1);
L = [];
Y0_current = Y0;

fprintf('=== PHUONG PHAP LUY THUA KET HOP XUONG THANG ===\n\n');
fprintf('Ma tran A ban dau:\n');
disp_matrix(A);
fprintf('Vector khoi tao Y0:\n');
disp_vector(Y0_current);
fprintf('Sai so epsilon = %.7f\n', epsilon);
fprintf('So buoc lap toi da N0 = %d\n\n', N0);

break_outer = false;

for i = 1:n
    fprintf('---------------------------------------------------\n');
    fprintf('VONG LAP LON i = %d (Tim tri rieng thu %d)\n', i, i);
    fprintf('---------------------------------------------------\n');
    fprintf('Ma tran A o dau vong lap %d:\n', i);
    disp_matrix(A);

    B_hist = zeros(n, N0 + 1);
    Y_hist = zeros(n, N0 + 1);
    B_hist(:, 1) = Y0_current;

    converged_type = 0;
    m_final = N0;

    for m = 1:N0
        Y_hist(:, m+1) = A * B_hist(:, m);

        max_abs_Y = max(abs(Y_hist(:, m+1)));
        if max_abs_Y == 0
            B_hist(:, m+1) = Y_hist(:, m+1);
        else
            B_hist(:, m+1) = Y_hist(:, m+1) / max_abs_Y;
        end

        fprintf('  Buoc m=%2d | Y^(%d)^T: [', m, m+1);
        for d = 1:n
            fprintf('%.7f', Y_hist(d, m+1));
            if d < n, fprintf(', '); end
        end
        fprintf('] | B^(%d)^T: [', m+1);
        for d = 1:n
            fprintf('%.7f', B_hist(d, m+1));
            if d < n, fprintf(', '); end
        end
        fprintf(']\n');

        if m >= 3
            diff1 = norm(B_hist(:, m) - B_hist(:, m-1), 'inf');
            diff2 = norm(B_hist(:, m) - B_hist(:, m-2), 'inf');

            if diff1 <= epsilon
                converged_type = 1;
                m_final = m;
                fprintf('  -> Hoi tu o buoc m = %d: Tri rieng thuc troi (sai lech = %.7f)\n', m, diff1);
                break;
            elseif diff2 <= epsilon
                converged_type = 2;
                m_final = m;
                fprintf('  -> Hoi tu o buoc m = %d: Hai tri rieng doi nhau (sai lech = %.7f)\n', m, diff2);
                break;
            end
        end

        if m == N0
            converged_type = 3;
            m_final = N0;
            fprintf('  -> Dat gioi han lap m = N0 = %d: Chuyen sang xu ly cap tri rieng phuc\n', N0);
        end
    end

    if converged_type == 1
        V = B_hist(:, m_final);
        [~, k] = max(abs(V));
        AV = Y_hist(:, m_final + 1);
        lambda = AV(k) / V(k);

        L = [L; lambda];
        fprintf('\n  => KET QUA VONG %d: Tim thay tri rieng thuc troi lambda = %.7f\n', i, lambda);

        A = Xuong_Thang(A, V);
        fprintf('  Ma tran A sau khi xuong thang:\n');
        disp_matrix(A);

        Y0_current = rand(n, 1);
        fprintf('  Cap nhat Y0 ngau nhien cho vong sau^T: [');
        for d = 1:n
            fprintf('%.7f', Y0_current(d));
            if d < n, fprintf(', '); end
        end
        fprintf(']\n\n');

    elseif converged_type == 2
        M = B_hist(:, m_final);
        AM = Y_hist(:, m_final + 1);
        A2M = A * AM;

        [~, k] = max(abs(AM));
        lambda = sqrt(A2M(k) / M(k));

        L = [L; lambda];
        fprintf('\n  => KET QUA VONG %d: Tim thay tri rieng lambda = %.7f\n', i, lambda);
        fprintf('  (Luu y nghiem doi -%.7f se tu dong duoc tim thay o vong sau)\n', lambda);

        V = AM + lambda * M;
        max_abs_V = max(abs(V));
        if max_abs_V > 0
            V = V / max_abs_V;
        end

        A = Xuong_Thang(A, V);
        fprintf('  Ma tran A sau khi xuong thang:\n');
        disp_matrix(A);

        Y0_current = rand(n, 1);
        fprintf('\n');

    elseif converged_type == 3
        M = B_hist(:, m_final);
        AM = Y_hist(:, m_final + 1);
        A2M = A * AM;

        Mat_sys = [AM(1), M(1); AM(2), M(2)];
        Vec_sys = [-A2M(1); -A2M(2)];

        if abs(det(Mat_sys)) < 1e-12 && n >= 3
            Mat_sys = [AM(2), M(2); AM(3), M(3)];
            Vec_sys = [-A2M(2); -A2M(3)];
        end

        coeffs = Mat_sys \ Vec_sys;
        b = coeffs(1);
        c = coeffs(2);

        fprintf('\n  Giai phuong trinh dac trung bac hai voi b = %.7f, c = %.7f\n', b, c);
        r_roots = roots([1, b, c]);
        lambda1 = r_roots(1);
        lambda2 = r_roots(2);

        L = [L; lambda1; lambda2];

        fprintf('  => KET QUA VONG %d: Tim thay cap tri rieng phuc lien hop:\n', i);
        if imag(lambda1) >= 0
            fprintf('     lambda_1 = %.7f + %.7fi\n', real(lambda1), imag(lambda1));
        else
            fprintf('     lambda_1 = %.7f - %.7fi\n', real(lambda1), abs(imag(lambda1)));
        end

        if imag(lambda2) >= 0
            fprintf('     lambda_2 = %.7f + %.7fi\n', real(lambda2), imag(lambda2));
        else
            fprintf('     lambda_2 = %.7f - %.7fi\n', real(lambda2), abs(imag(lambda2)));
        end

        fprintf('  Ket thuc som thuat toan vi da tim thay du nghiem tu cap phuc.\n\n');
        break_outer = true;
    end

    if break_outer
        break;
    end
end

% =========================================================================
% 4. TONG HOP SO LIEU QUAN TRONG CUOI CUNG
% =========================================================================
fprintf('===================================================\n');
fprintf('TONG HOP SO LIEU QUAN TRONG CUOI CUNG\n');
fprintf('===================================================\n');
fprintf('Ma tran ban dau A_init:\n');
disp_matrix(A_init);
fprintf('Tong so tri rieng da tim thay: %d\n', length(L));
fprintf('Danh sach cac tri rieng chi tiet:\n');
for idx = 1:length(L)
    val = L(idx);
    if imag(val) == 0
        fprintf('  Tri rieng %d: %.7f\n', idx, real(val));
    elseif imag(val) >= 0
        fprintf('  Tri rieng %d: %.7f + %.7fi\n', idx, real(val), imag(val));
    else
        fprintf('  Tri rieng %d: %.7f - %.7fi\n', idx, real(val), abs(imag(val)));
    end
end
fprintf('===================================================\n');
