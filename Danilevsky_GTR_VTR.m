clc;
clear;

% =========================================================================
% CAU HINH BAI TOAN
% =========================================================================
% Ma tran input vuong A
A = [1, 2, 3, 4;
     2, 1, 4, 3;
     3, 4, 1, 2;
     4, 3, 2, 1];

% =========================================================================
% CAC HAM HO TRO IN AN (TIENG VIET KHONG DAU, 7 CHU SO THAP PHAN)
% =========================================================================
function in_ma_tran(M, ten)
    fprintf('%s =\n', ten);
    [r, c] = size(M);
    for i = 1:r
        for j = 1:c
            re = real(M(i,j));
            im = imag(M(i,j));
            if abs(im) > 1e-9
                if im > 0
                    fprintf('%12.7f + %12.7fi ', re, im);
                else
                    fprintf('%12.7f - %12.7fi ', re, abs(im));
                end
            else
                fprintf('%12.7f ', re);
            end
        end
        fprintf('\n');
    end
    fprintf('\n');
end

% =========================================================================
% CAC HAM BO TRO TOAN HOC
% =========================================================================
function [A_new, M] = findSimpleA(A, k)
    m = size(A, 1);
    M = eye(m);
    for j = 1:m
        if j ~= k-1
            M(k-1, j) = -A(k, j) / A(k, k-1);
        else
            M(k-1, j) = 1 / A(k, k-1);
        end
    end
    A_new = inv(M) * A * M;
end

function [L, V] = solve_frobenius(B, S_part)
    m = size(B, 1);
    p_coeffs = zeros(1, m + 1);
    p_coeffs(1) = 1;
    for i = 1:m
        p_coeffs(i+1) = -B(1, i);
    end

    L = roots(p_coeffs);
    V = zeros(size(S_part, 1), length(L));

    for i = 1:length(L)
        lam = L(i);
        u = zeros(m, 1);
        for idx = 1:m
            u(idx) = lam^(m - idx);
        end
        v = S_part * u;
        v = v / norm(v, 2);
        V(:, i) = v;
    end
end

% =========================================================================
% THUAT TOAN CHINH: Danilevsky_Main
% =========================================================================
fprintf('--- BAT DAU THUAT TOAN ---\n\n');
in_ma_tran(A, 'Ma tran A ban dau');

n = size(A, 1);
m = n;
k = n;
S = eye(n);
L_all = [];
V_all = [];

buoc = 1;

while k > 1
    fprintf('--- Buoc lap %d, k = %d, m = %d ---\n', buoc, k, m);

    if abs(A(k, k-1)) > 1e-9
        fprintf('TH 1: Khu binh thuong\n');
        [A, M_k] = findSimpleA(A, k);
        S(:, 1:m) = S(:, 1:m) * M_k;

        in_ma_tran(M_k, 'Ma tran bien doi M_k');
        in_ma_tran(A, 'Ma tran A cap nhat');
        in_ma_tran(S, 'Ma tran toan cuc S cap nhat');
        k = k - 1;

    else
        j_doi = 0;
        for j = 1:(k-2)
            if abs(A(k, j)) > 1e-9
                j_doi = j;
                break;
            end
        end

        if j_doi > 0
            fprintf('TH 2: Doi cho cot %d va %d\n', j_doi, k-1);
            % Doi hang va cot tren A
            A([j_doi, k-1], :) = A([k-1, j_doi], :);
            A(:, [j_doi, k-1]) = A(:, [k-1, j_doi]);

            % Doi cot tren S
            S(:, [j_doi, k-1]) = S(:, [k-1, j_doi]);

            in_ma_tran(A, 'Ma tran A sau hoan vi');
            in_ma_tran(S, 'Ma tran S sau hoan vi');

        else
            fprintf('TH 3: Tach khoi\n');

            % 1. Chuan hoa khoi duoi
            for j = k:(m-2)
                M_tmp = eye(m);
                M_tmp(1:k, j+1) = -A(1:k, j);
                invM = inv(M_tmp);
                A = M_tmp * A * invM;
                S(:, 1:m) = S(:, 1:m) * invM;
            end
            if k <= m-2
                fprintf('Da chuan hoa khoi duoi.\n');
            end

            % 2. Kiem tra lien ket
            dich_hang = false;
            for j = (k-1):-1:1
                if abs(A(j, m)) > 1e-9
                    fprintf('Phat hien lien ket tai hang %d, day xuong day khoi tren.\n', j);
                    P = eye(m);
                    P([k-1, j], :) = P([j, k-1], :);

                    A = P * A * P';
                    S(:, 1:m) = S(:, 1:m) * P';
                    k = m;
                    dich_hang = true;

                    in_ma_tran(A, 'Ma tran A cap nhat sau khi day hang');
                    break;
                end
            end

            % 3. Giai va thu hep khoi
            if ~dich_hang
                fprintf('Tien hanh giai khoi Frobenius va thu hep.\n');
                B = A(k:m, k:m);
                in_ma_tran(B, 'Khoi Frobenius B');

                [L_part, V_part] = solve_frobenius(B, S(:, k:m));
                L_all = [L_all; L_part];
                V_all = [V_all, V_part];

                m = k - 1;
                A = A(1:m, 1:m);
                k = m;

                in_ma_tran(A, 'Ma tran A thu hep');
            end
        end
    end
    buoc = buoc + 1;
end

% =========================================================================
% BUOC 3: XU LY KHOI CUOI CUNG
% =========================================================================
fprintf('\n--- XU LY KHOI CUOI CUNG ---\n');
in_ma_tran(A, 'Khoi Frobenius cuoi cung');
[L_fin, V_fin] = solve_frobenius(A, S(:, 1:m));

L_all = [L_all; L_fin];
V_all = [V_all, V_fin];

% =========================================================================
% BUOC 4: TONG HOP KET QUA
% =========================================================================
fprintf('\n=========================================================================\n');
fprintf('KET QUA CUOI CUNG\n');
fprintf('=========================================================================\n');
in_ma_tran(L_all, 'Danh sach tri rieng list_ev');
in_ma_tran(V_all, 'Danh sach vector rieng tuong ung list_vec (moi cot la 1 vector)');
