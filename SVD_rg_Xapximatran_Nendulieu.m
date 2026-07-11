% =========================================================================
% CHƯƠNG TRÌNH PHÂN TÍCH SVD RÚT GỌN VÀ XẤP XỈ MA TRẬN HẠNG NHỎ
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
p = 1;             % Số lượng giá trị kỳ dị giữ lại để xấp xỉ (hạng mục tiêu)

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

    x = ones(n, 1);
    x = x / norm(x);

    w = x;
    thoat_hoan_toan = false;

    for j = 1:N0
        Y = M * x;

        if norm(Y) <= epsilon
            r = i - 1;
            thoat_hoan_toan = true;
            fprintf('  -> ||Y|| <= epsilon. Hang thuc te r = %d. Thoat hoan toan vong lap %d.\n', r, i);
            break;
        end

        x_new = Y / norm(Y);

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
    fprintf('\n');

    U_res = [U_res, u];
    S_res = [S_res; sigma];
    V_res = [V_res, v];
end

% B3: Tổng hợp kết quả
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

A_approx = U * Sigma * V';
sai_so = norm(A - A_approx, 'fro');
fprintf('\nKiem tra sai so ||A - U*Sigma*V^T||_F = %.7f\n', sai_so);

% =========================================================================
% BỔ SUNG: XẤP XỈ MA TRẬN HẠNG NHỎ & NÉN DỮ LIỆU
% =========================================================================
fprintf('\n=== GIAI BAI TOAN XAP XI MA TRAN (NEN DU LIEU) ===\n');

% Xử lý trường hợp p lớn hơn hạng thực tế r
if p > r
    p = r;
end
fprintf('Chon muc do nen (hang p) = %d\n\n', p);

% 1. Xây dựng ma trận xấp xỉ A_hat
A_hat = zeros(m, n);
for k = 1:p
    A_hat = A_hat + S_res(k) * U_res(:, k) * V_res(:, k)';
end
fprintf('Ma tran xap xi A_hat (giai nen de hien thi):\n');
in_ma_tran(A_hat);

% 2. Đánh giá lượng thông tin theo chuẩn Frobenius
tong_bp_sigma_goc = sum(S_res.^2);
tong_bp_sigma_giu = sum(S_res(1:p).^2);
tong_bp_sigma_bo  = 0;

if p < r
    tong_bp_sigma_bo = sum(S_res(p+1:r).^2);
end

phan_tram_giu_lai = sqrt(tong_bp_sigma_giu / tong_bp_sigma_goc) * 100;
phan_tram_bo_di   = sqrt(tong_bp_sigma_bo / tong_bp_sigma_goc) * 100;

fprintf('\n--- Danh gia luong thong tin ---\n');
fprintf('Phan tram thong tin giu lai: %.7f %%\n', phan_tram_giu_lai);
fprintf('Phan tram thong tin bo di (sai so tuong doi): %.7f %%\n', phan_tram_bo_di);

% 3. Sai số tuyệt đối (Chuẩn phổ)
if p < r
    sai_so_tuyet_doi = S_res(p+1);
else
    sai_so_tuyet_doi = 0;
end
fprintf('Sai so xap xi tuyet doi ||A - A_hat||_2 = %.7f\n', sai_so_tuyet_doi);

% 4. Đánh giá hiệu quả nén lưu trữ
dung_luong_goc = m * n;
dung_luong_nen = p * (1 + m + n);

fprintf('\n--- Hieu qua nen luu tru ---\n');
fprintf('So luong phan tu can luu (goc): %d\n', dung_luong_goc);
fprintf('So luong phan tu can luu (nen): %d\n', dung_luong_nen);
fprintf('Ty le dung luong (Nen / Goc): %.7f\n', dung_luong_nen / dung_luong_goc);
