clc;
clear;

% =========================================================================
% CAU HINH BAI TOAN
% =========================================================================
% Ma tran input vuong A (Ban co the thay doi ma tran tai day)
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
            fprintf('%12.7f ', M(i,j));
        end
        fprintf('\n');
    end
    fprintf('\n');
end

function in_da_thuc(p, ten)
    fprintf('%s = ', ten);
    n = length(p);
    for i = 1:n
        bac = n - i;
        if i > 1 && p(i) >= 0
            fprintf('+ ');
        elseif p(i) < 0
            fprintf('- ');
        end

        val = abs(p(i));
        if bac > 1
            fprintf('%.7f*L^%d ', val, bac);
        elseif bac == 1
            fprintf('%.7f*L ', val);
        else
            fprintf('%.7f\n', val);
        end
    end
    fprintf('\n');
end

% =========================================================================
% THUAT TOAN DANILEVSKY
% =========================================================================
fprintf('--- BAT DAU THUAT TOAN ---\n\n');
in_ma_tran(A, 'Ma tran A ban dau');

n = size(A, 1);
m = n;
k = n;
P_A = [1]; % Khoi tao da thuc P_A(L) = 1 (L la lambda)

buoc = 1;

while k > 1
    fprintf('--- Buoc lap %d (k = %d, m = %d) ---\n', buoc, k, m);

    if abs(A(k, k-1)) > 1e-9 % TH 1: Khu binh thuong
        fprintf('TH 1: Khu binh thuong (A(k, k-1) != 0)\n');
        M = eye(m);
        for j = 1:m
            if j ~= k-1
                M(k-1, j) = -A(k, j) / A(k, k-1);
            else
                M(k-1, j) = 1 / A(k, k-1);
            end
        end
        in_ma_tran(M, 'Ma tran M');

        invM = inv(M);
        in_ma_tran(invM, 'Ma tran nghich dao M^-1');

        A = invM * A * M;
        in_ma_tran(A, 'Ma tran A sau khi cap nhat');

        k = k - 1;

    else % Kiem tra xem co the roi vao TH 2 hay TH 3
        j_doi = 0;
        for j = 1:(k-2)
            if abs(A(k, j)) > 1e-9
                j_doi = j;
                break;
            end
        end

        if j_doi > 0 % TH 2: Doi cho cot
            fprintf('TH 2: Doi cho cot (A(k, k-1) = 0 va A(k, %d) != 0)\n', j_doi);

            % Doi cho hang j va k-1
            temp_row = A(j_doi, :);
            A(j_doi, :) = A(k-1, :);
            A(k-1, :) = temp_row;

            % Doi cho cot j va k-1
            temp_col = A(:, j_doi);
            A(:, j_doi) = A(:, k-1);
            A(:, k-1) = temp_col;

            in_ma_tran(A, 'Ma tran A sau khi doi cho');
            % k giu nguyen theo thuat toan

        else % TH 3: Tach khoi
            fprintf('TH 3: Tach khoi (Moi A(k, j) = 0 voi j <= k-2)\n');

            % 1. Chuan hoa khoi duoi
            for j = k:(m-2)
                M = eye(m);
                M(1:k, j+1) = -A(1:k, j);
                A = M * A * inv(M);
            end
            if k <= m-2
                in_ma_tran(A, 'Ma tran A sau khi chuan hoa khoi duoi');
            end

            % 2. Kiem tra lien ket
            dich_hang = false;
            for j = (k-1):-1:1
                if abs(A(j, m)) > 1e-9
                    fprintf('Phat hien lien ket tai A(%d, %d) != 0. Thuc hien day hang.\n', j, m);
                    P = eye(m);
                    % Day hang k-1 xuong hoan vi voi j
                    temp = P(k-1, :);
                    P(k-1, :) = P(j, :);
                    P(j, :) = temp;

                    A = P * A * P';
                    k = m;
                    dich_hang = true;
                    in_ma_tran(A, 'Ma tran A sau hoan vi P');
                    break;
                end
            end

            % 3. Thu hep khoi (neu khong dich hang)
            if ~dich_hang
                fprintf('Khong co lien ket, thuc hien thu hep khoi.\n');
                r = m - k + 1;
                B = A(k:m, k:m);
                in_ma_tran(B, 'Khoi Frobenius B trich xuat');

                % Da thuc p_B = L^r - B(1,1)L^(r-1) - ... - B(1,r)
                p_B = zeros(1, r + 1);
                p_B(1) = 1;
                for i = 1:r
                    p_B(i+1) = -B(1, i);
                end
                in_da_thuc(p_B, 'Da thuc dac trung khoi B (p_B)');

                % Tich luy da thuc bang phep nhan chap (convolution)
                P_A = conv(P_A, p_B);
                in_da_thuc(P_A, 'Da thuc P_A tich luy');

                % Thu hep ma tran
                m = k - 1;
                A = A(1:m, 1:m);
                k = m;
                in_ma_tran(A, 'Ma tran A sau khi thu hep');
            end
        end
    end
    buoc = buoc + 1;
end

% =========================================================================
% BUOC 3: XU LY KHOI CUOI CUNG
% =========================================================================
fprintf('\n--- XU LY KHOI CUOI CUNG ---\n');
F = A(1:m, 1:m);
in_ma_tran(F, 'Khoi Frobenius cuoi cung F');

p_F = zeros(1, m + 1);
p_F(1) = 1;
for i = 1:m
    p_F(i+1) = -F(1, i);
end
in_da_thuc(p_F, 'Da thuc dac trung khoi F (p_F)');

P_A = conv(P_A, p_F);

% =========================================================================
% BUOC 4: TONG HOP & IN KET QUA
% =========================================================================
fprintf('\n=========================================================================\n');
fprintf('TONG HOP KET QUA QUAN TRONG\n');
fprintf('=========================================================================\n');
in_da_thuc(P_A, 'DA THUC DAC TRUNG CUA MA TRAN A (P_A(L))');
