clc;
format short;
% ==========================================
% CAU HINH THAM SO DAU VAO
% ==========================================
% Ma tran A
A = [ 5.1 4.1 3.1 2.1;
      4.1 6.1 2.1 5.1;
      3.1 2.1 5.1 4.1;
      2.1 5.1 4.1 6.1];

% Tu dong lay kich thuoc ma tran (khong can nhap thu cong n)
n = size(A, 1);

E = 0.00000001;          % Sai so mong muon
max_iter = 1000;         % So buoc lap toi da

% Vector khoi tao Y
Y = [0;
     1;
     0;
     0];

% ==========================================
% CAC HAM PHU TRO (FUNCTIONS)
% ==========================================

% Chuan hoa vector
function matrix = Chuan_Hoa(M)
    [~, idx] = max(abs(M));
    c = M(idx);
    matrix = M / c;
end

% Kiem tra sai so giua 2 vector
function err = Kiem_tra(B, m, n_idx)
    A_diff = B{m} - B{n_idx};
    err = max(abs(A_diff));
end

% Xu ly (Vong lap luy thua tim vector hoi tu)
function [m, TH, B] = Xu_li(A, Y, n, E, max_iter)
    B = cell(1, 1);
    B{1} = Y;
    Z = zeros(n, 1);

    B{2} = A * B{1};
    B{3} = A * B{2};
    B{4} = A * B{3};

    m = 4;
    TH = 4;

    while true
        F = round(B{end} * 10000) / 10000;

        if all(F == Z)
            error('Vector hoi tu ve 0, thuat toan dung.');
        end

        % Dieu kien dung: Dat den so vong lap toi da
        if m >= max_iter
            break;
        end

        % Dieu kien dung: Dat duoc sai so mong muon
        if Kiem_tra(B, m - 1, m - 2) <= E
            TH = 1;
            break;
        end

        if Kiem_tra(B, m - 1, m - 3) <= E
            TH = 3;
            break;
        end

        m = m + 1;
        B{m} = Chuan_Hoa(A * B{m - 1});
    end
end

% Chuyen ma tran (Phuong phap xuong thang)
function A_new = Chuyen_MT(A, VTR, n)
    V = VTR(:);
    M = eye(n);

    maxi_1 = max(V);
    index = find(V == maxi_1, 1);

    E_col = M(:, index) - V;
    E_col(index) = 0;
    M(:, index) = E_col;

    A_new = M * A;
end

% ==========================================
% CHUONG TRINH CHINH
% ==========================================
fprintf('Ma tran co: %d x %d\n', n, n);
disp('Ma tran A la: ');
disp(A);
disp(repmat('_', 1, 100));

for i = 1:n

    [m_iters, TH, B] = Xu_li(A, Y, n, E, max_iter);

    if TH == 1

        fprintf('=> Tim duoc lamda theo TRUONG HOP %d\n\n', TH);

        VTR = B{end};
        VTR_T = VTR';
        AX = A * VTR;

        % Tìm vị trí phần tử có trị tuyệt đối lớn nhất trong VTR
        [~, idx] = max(abs(VTR));

        % Lấy thành phần tương ứng của AX chia cho VTR
        lamda = AX(idx) / VTR(idx);
        lamda = round(lamda * 10000) / 10000;

        fprintf('Sau %d lan lap\n', m_iters - 1);
        fprintf('lamda_%d la: %f\n', i, lamda);
        fprintf('Vector rieng ung voi lamda_%d la: \n', i);
        disp(VTR);

        A = Chuyen_MT(A, VTR, n);

        fprintf('\nMa tran sau khi xuong thang:\n');
        disp(A);

        Y = rand(n, 1);
        disp(repmat('_', 1, 100));

    elseif TH == 3

        fprintf('=> Tim duoc lamda theo TRUONG HOP %d\n\n', TH);

        AmY = A * B{end};
        Am1Y = A * AmY;
        Am2Y = A * Am1Y;

        lamda_N = max(Am2Y) / max(AmY);
        lamda_N = sqrt(lamda_N);

        VTR_N = Am1Y + lamda_N * AmY;
        VTR_N = Chuan_Hoa(VTR_N);

        fprintf('Sau %d lan lap\n', m_iters - 1);
        fprintf('lamda_%d la: %f\n', i, round(lamda_N * 10000) / 10000);
        fprintf('Vector rieng ung voi lamda_%d la: \n', i);
        disp(VTR_N);

        A = Chuyen_MT(A, VTR_N, n);

        fprintf('\nMa tran sau khi xuong thang:\n');
        disp(A);

        Y = rand(n, 1);
        disp(repmat('_', 1, 100));

    elseif TH == 4

        fprintf('=> Tim duoc lamda theo TRUONG HOP %d\n\n', TH);

        M = B{end};
        AmY = A * M;
        Am1Y = A * AmY;

        a_1 = Am1Y(1); a_2 = Am1Y(2);
        b_1 = AmY(1);  b_2 = AmY(2);
        c_1 = M(1);    c_2 = M(2);

        a_coef = 1;
        b_coef = (a_1 * c_2 - c_1 * a_2) / (c_1 * b_2 - b_1 * c_2);
        c_coef = (b_1 * a_2 - a_1 * b_2) / (c_1 * b_2 - b_1 * c_2);

        fprintf('Phuong trinh tim duoc: Z^2 + %fZ + %f = 0 \n\n', b_coef, c_coef);

        denta = b_coef^2 - 4 * a_coef * c_coef;

        if denta >= 0
            disp('Phuong trinh co denta >= 0');
        else
            lamda_1 = complex(-b_coef / (2 * a_coef), ...
                              -sqrt(abs(denta)) / (2 * a_coef));

            lamda_2 = complex(-b_coef / (2 * a_coef), ...
                               sqrt(abs(denta)) / (2 * a_coef));

            VTR_1 = Am1Y - lamda_1 * AmY;
            VTR_2 = Am1Y - lamda_2 * AmY;

            fprintf('lamda_%d la:\n', i);
            disp(lamda_1);

            fprintf('Vector rieng ung voi lamda_%d la: \n', i);
            disp(VTR_1);

            disp(repmat('*', 1, 100));

            fprintf('lamda_%d la:\n', i + 1);
            disp(lamda_2);

            fprintf('Vector rieng ung voi lamda_%d la: \n', i + 1);
            disp(VTR_2);
        end

        break;
    end
end
