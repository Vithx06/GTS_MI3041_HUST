% =========================================================================
% CHƯƠNG TRÌNH PHÂN TÍCH SVD RÚT GỌN & MA TRẬN GIẢ NGHỊCH ĐẢO
% =========================================================================

clc; clear; close all;

% --- CẤU HÌNH BÀI TOÁN ---
% Ma trận đầu vào A (Bạn có thể thay đổi ma trận tại đây)
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
r = min(m, n);

U_res = []; % Lưu trữ các vector u
S_res = []; % Lưu trữ các giá trị sigma
V_res = []; % Lưu trữ các vector v

fprintf('Khoi tao M(1) = A^T * A:\n'); in_ma_tran(M);
fprintf('Hang toi da r = %d\n\n', r);

% B2: Vòng lặp chính
for i = 1:min(m, n)
    fprintf('=== BUOC %d ===\n', i);

    % Có thể dùng vector ngẫu nhiên hoặc vector toàn 1 rồi chuẩn hóa.
    x = ones(n, 1);
    x = x / norm(x);

    % 2.2: Lặp lũy thừa tìm vector riêng
    w = x;
    thoat_hoan_toan = false;

    for j = 1:N0
        Y = M * x;

        % Kiểm tra kết thúc sớm
        if norm(Y) <= epsilon
            r = i - 1;
            thoat_hoan_toan = true;
            fprintf('  -> ||Y|| <= epsilon. Hang thuc te r = %d. Thoat hoan toan vong lap %d.\n', r, i);
            break;
        end

        % Chuẩn hóa
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

    % 2.3: Tính trị riêng và giá trị kì dị
    lambda = w' * M * w;
    sigma = sqrt(lambda);
    fprintf('  Tri rieng lambda_%d = %.7f\n', i, lambda);
    fprintf('  Gia tri ki di sigma_%d = %.7f\n', i, sigma);

    % 2.4: Vector kì dị phải
    v = w;
    fprintf('  Vector ki di phai v_%d =\n', i);
    in_ma_tran(v);

    % 2.5: Vector kì dị trái
    u = (A * v) / sigma;
    fprintf('  Vector ki di trai u_%d =\n', i);
    in_ma_tran(u);

    % 2.6: Lập ma trận xuống thang
    M = M - lambda * (v * v');
    fprintf('  Ma tran xuong thang M^(%d) =\n', i + 1);
    in_ma_tran(M);
    fprintf('\n');

    % Lưu kết quả
    U_res = [U_res, u];
    S_res = [S_res; sigma];
    V_res = [V_res, v];
end

% B3: Tổng hợp kết quả SVD rút gọn
fprintf('=== TONG HOP KET QUA (A = U * Sigma * V^T) ===\n');

U = U_res;
Sigma = diag(S_res);
V = V_res;

fprintf('Ma tran U (%d x %d):\n', size(U, 1), size(U, 2));
in_ma_tran(U);

fprintf('\nMa tran Sigma (%d x %d):\n', size(Sigma, 1), size(Sigma, 2));
in_ma_tran(Sigma);

fprintf('\nMa tran V (%d x %d):\n', size(V, 1), size(V, 2));
in_ma_tran(V);

% Kiểm tra lại sai số phục hồi A
A_approx = U * Sigma * V';
sai_so = norm(A - A_approx, 'fro');
fprintf('\nKiem tra sai so ||A - U*Sigma*V^T||_F = %.7f\n', sai_so);

% =========================================================================
% BỔ SUNG: TÍNH MA TRẬN GIẢ NGHỊCH ĐẢO (MOORE-PENROSE PSEUDOINVERSE)
% =========================================================================
fprintf('\n=== TINH MA TRAN GIA NGHICH DAO A^dagger ===\n');

% Tạo ma trận Sigma_dagger bằng cách lấy nghịch đảo các giá trị kì dị
S_dagger_res = zeros(size(S_res));
for k = 1:length(S_res)
    if S_res(k) > epsilon % Đảm bảo không chia cho 0 với các giá trị cực nhỏ
        S_dagger_res(k) = 1 / S_res(k);
    end
end

Sigma_dagger = diag(S_dagger_res);

fprintf('Ma tran Sigma^dagger (%d x %d):\n', size(Sigma_dagger, 1), size(Sigma_dagger, 2));
in_ma_tran(Sigma_dagger);

% Tính A_dagger = V * Sigma_dagger * U^T
A_dagger = V * Sigma_dagger * U';

fprintf('\nMa tran gia nghich dao A^dagger (%d x %d):\n', size(A_dagger, 1), size(A_dagger, 2));
in_ma_tran(A_dagger);

% Có thể so sánh nhanh với hàm giả nghịch đảo có sẵn của Octave: pinv(A)
% fprintf('\nKiem tra so sanh voi ham pinv(A) cua Octave:\n');
% in_ma_tran(pinv(A));
