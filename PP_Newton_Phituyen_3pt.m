% CHƯƠNG TRÌNH PHƯƠNG PHÁP TIẾP TUYẾN (NEWTON) CHO HỆ 3 PHƯƠNG TRÌNH
clc; clear;

% B1: Nhập các hàm số f, g, h, vector xấp xỉ đầu và số bước lặp N
% Hệ phương trình: f(x,y,z)=0, g(x,y,z)=0, h(x,y,z)=0
f = @(x, y, z) x^2 + y^2 + z^2 - 9;
g = @(x, y, z) x*y*z - 1;
h = @(x, y, z) x + y - z^2;

X0 = [2.5; 1.5; 1.0];   % Vector cột xấp xỉ đầu (x0, y0, z0)^T
N = 5;                  % Số bước lặp yêu cầu

% B2: Thiết lập Ma trận Jacobi J và Vector hàm F
% Tính các đạo hàm riêng (9 đạo hàm cho hệ 3 ẩn)
df_dx = @(x, y, z) 2*x;   df_dy = @(x, y, z) 2*y;   df_dz = @(x, y, z) 2*z;
dg_dx = @(x, y, z) y*z;   dg_dy = @(x, y, z) x*z;   dg_dz = @(x, y, z) x*y;
dh_dx = @(x, y, z) 1;     dh_dy = @(x, y, z) 1;     dh_dz = @(x, y, z) -2*z;

% Vector hàm F(X)
F = @(X) [f(X(1), X(2), X(3));
          g(X(1), X(2), X(3));
          h(X(1), X(2), X(3))];

% Ma trận Jacobi J(X) kích thước 3x3
J = @(X) [df_dx(X(1), X(2), X(3)), df_dy(X(1), X(2), X(3)), df_dz(X(1), X(2), X(3));
          dg_dx(X(1), X(2), X(3)), dg_dy(X(1), X(2), X(3)), dg_dz(X(1), X(2), X(3));
          dh_dx(X(1), X(2), X(3)), dh_dy(X(1), X(2), X(3)), dh_dz(X(1), X(2), X(3))];

% B3: Vòng lặp i chạy từ 0 đến N - 1
X_old = X0;
fprintf('Bat dau qua trinh lap Newton cho he 3 phuong trinh:\n');
fprintf('Buoc 0: x = %.10f, y = %.10f, z = %.10f\n', X0(1), X0(2), X0(3));

for i = 0:(N-1)
    % Tính giá trị Jacobi và Vector hàm tại điểm hiện tại
    Ji = J(X_old);
    Fi = F(X_old);

    % In ma trận Jacobi để kiểm tra (tùy chọn)
    % fprintf('\nMa tran Jacobi buoc %d:\n', i); disp(Ji);

    % 3.1: Kiểm tra điều kiện định thức det(J) != 0
    if det(Ji) == 0
        error('Thong bao: Ma tran Jacobi suy bien (det = 0). Ket thuc thuat toan.');
    end

    % 3.2: Tính X_{i+1} = X_i - J^-1 * F(X_i)
    % Sử dụng toán tử \ thay cho inv() để tối ưu hóa tính toán
    X_new = X_old - Ji \ Fi;

    % Cập nhật cho bước sau
    X_old = X_new;

    % In kết quả mỗi bước lặp
    fprintf('Buoc %d: x = %.10f, y = %.10f, z = %.10f\n', i + 1, X_new(1), X_new(2), X_new(3));
end

% B4: In ra kết quả cuối cùng
fprintf('\n--- KET QUA CUOI CUNG ---\n');
fprintf('Vec-to nghiem gan dung X = [%.10f; %.10f; %.10f]\n', X_old(1), X_old(2), X_old(3));
fprintf('Gia tri cac ham tai nghiem:\n f = %.2e\n g = %.2e\n h = %.2e\n', ...
        f(X_old(1), X_old(2), X_old(3)), ...
        g(X_old(1), X_old(2), X_old(3)), ...
        h(X_old(1), X_old(2), X_old(3)));
