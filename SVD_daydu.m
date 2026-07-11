% =========================================================================
% CHƯƠNG TRÌNH PHÂN TÍCH SVD ĐẦY ĐỦ (PHƯƠNG PHÁP LŨY THỪA & XUỐNG THANG)
% =========================================================================

clc; clear; close all;

% --- CẤU HÌNH BÀI TOÁN ---
A = [
    3   3   7;
    3   9  10;
    8   2  10;
    2  10   9
];

epsilon = 1e-6;    % Sai số dừng
N0 = 1000;         % Số lần lặp tối đa

% Hàm hỗ trợ in ma trận/vector với 7 chữ số sau dấu phẩy
function in_ma_tran(mat)
    [rows, cols] = size(mat);
    for i = 1:rows
        for j = 1:cols
            fprintf('%12.7f ', mat(i, j));
        end
        fprintf('\n');
    end
end

% Hàm hỗ trợ hoàn thiện cơ sở trực chuẩn bằng Gram-Schmidt
function Base_out = Gram_Schmidt_hoan_thien(Base_in, dim, eps_tol)
    Base_out = Base_in;
    num_vec = size(Base_out, 2);

    if num_vec == dim
        return; % Đã đủ chiều, không cần thêm
    end

    I = eye(dim); % Lấy các vector cơ sở chuẩn làm ứng viên
    for k = 1:dim
        e = I(:, k);
        % Trực giao hóa e đối với tất cả các vector hiện có trong Base_out
        for j = 1:size(Base_out, 2)
            e = e - (Base_out(:, j)' * e) * Base_out(:, j);
        end

        % Nếu vector dư có độ dài đáng kể, chuẩn hóa và thêm vào cơ sở
        if norm(e) > eps_tol
            e = e / norm(e);
            Base_out = [Base_out, e];
            % Dừng lại nếu đã đủ số chiều
            if size(Base_out, 2) == dim
                break;
            end
        end
    end
end

% =========================================================================
% BẮT ĐẦU THUẬT TOÁN
% =========================================================================

fprintf('--- DỮ LIỆU ĐẦU VÀO ---\n');
fprintf('Ma tran A:\n'); in_ma_tran(A);
fprintf('Sai so dung epsilon = %.7f\n', epsilon);
fprintf('So lan lap toi da N0 = %d\n\n', N0);

% B1: Khởi tạo
[m, n] = size(A);
M = A' * A;
r_max = min(m, n);
r = r_max;

U_res = [];
S_res = [];
V_res = [];

fprintf('Khoi tao M(1) = A^T * A:\n'); in_ma_tran(M);
fprintf('Hang toi da r_max = %d\n\n', r_max);

% B2: Vòng lặp chính tìm các thành phần khác 0
for i = 1:r_max
    fprintf('=== BUOC %d (Tim thanh phan khac 0) ===\n', i);

    x = ones(n, 1);
    x = x / norm(x);

    w = x;
    thoat_hoan_toan = false;

    for j = 1:N0
        Y = M * x;

        % Kiểm tra kết thúc sớm (Hạng thực tế nhỏ hơn min(m,n))
        if norm(Y) <= epsilon
            r = i - 1;
            thoat_hoan_toan = true;
            fprintf('  -> ||Y|| <= epsilon. Hang thuc te r = %d. Thoat hoan toan vong lap %d.\n', r, i);
            break;
        end

        x_new = Y / norm(Y);

        % Kiểm tra hội tụ
        if norm(x_new - x) <= epsilon
            w = x_new;
            fprintf('  -> Hoi tu tai lan lap j = %d\n', j);
            break;
        end
        x = x_new;
    end

    if thoat_hoan_toan
        break;
    end

    lambda = w' * M * w;
    sigma = sqrt(lambda);
    fprintf('  Tri rieng lambda_%d = %.7f\n', i, lambda);
    fprintf('  Gia tri ki di sigma_%d = %.7f\n', i, sigma);

    v = w;
    fprintf('  Vector ki di phai v_%d =\n', i);
    in_ma_tran(v);

    u = (A * v) / sigma;
    fprintf('  Vector ki di trai u_%d =\n', i);
    in_ma_tran(u);

    M = M - lambda * (v * v');
    fprintf('  Ma tran xuong thang M^(%d) =\n', i + 1);
    in_ma_tran(M);

    U_res = [U_res, u];
    S_res = [S_res; sigma];
    V_res = [V_res, v];
    fprintf('\n');
end

% B3: Xây dựng ma trận kì dị phải V (n x n)
fprintf('=== BUOC 3: Xay dung ma tran ki di phai V (%d x %d) ===\n', n, n);
V = Gram_Schmidt_hoan_thien(V_res, n, epsilon);
in_ma_tran(V);
fprintf('\n');

% B4: Xây dựng ma trận kì dị trái U (m x m)
fprintf('=== BUOC 4: Xay dung ma tran ki di trai U (%d x %d) ===\n', m, m);
U = Gram_Schmidt_hoan_thien(U_res, m, epsilon);
in_ma_tran(U);
fprintf('\n');

% B5: Xây dựng ma trận đường chéo Sigma (m x n)
fprintf('=== BUOC 5: Xay dung ma tran duong cheo Sigma (%d x %d) ===\n', m, n);
Sigma = zeros(m, n);
for k = 1:r
    Sigma(k, k) = S_res(k);
end
in_ma_tran(Sigma);
fprintf('\n');

% B6: Xuất kết quả phân tích SVD đầy đủ
fprintf('=== BUOC 6: TONG HOP KET QUA SVD DAY DU (A = U * Sigma * V^T) ===\n');
A_approx = U * Sigma * V';
sai_so = norm(A - A_approx, 'fro');
fprintf('Kiem tra sai so ||A - U*Sigma*V^T||_F = %.7f\n', sai_so);
