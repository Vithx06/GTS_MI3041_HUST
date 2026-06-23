% CHƯƠNG TRÌNH PHƯƠNG PHÁP TIẾP TUYẾN (NEWTON) CHO HỆ 2 PHƯƠNG TRÌNH
clc; clear;

% B1: Nhập hàm số, vector xấp xỉ đầu và số bước lặp N
% Hệ phương trình: f(x,y) = 0 và g(x,y) = 0
f = @(x, y) x^2+2*y^2-9;
g = @(x, y) 5*y^2+x*y-4;

X0 = [-1.1; 3];         % Vector cột xấp xỉ đầu (x0, y0)^T
N = 5;                % Số bước lặp yêu cầu

% B2: Thiết lập Vector hàm F và Ma trận Jacobi J
% Nhập các đạo hàm riêng (Tính tay hoặc dùng công thức)
df_dx = @(x, y) 2*x;
df_dy = @(x, y) 4*y;
dg_dx = @(x, y) 10*x+y;
dg_dy = @(x, y) x;

% Vector hàm F(X)
F = @(X) [f(X(1), X(2));
          g(X(1), X(2))];

% Ma trận Jacobi J(X)
J = @(X) [df_dx(X(1), X(2)), df_dy(X(1), X(2));
          dg_dx(X(1), X(2)), dg_dy(X(1), X(2))];

% B3: Vòng lặp i chạy từ 0 đến N - 1
X_old = X0;
fprintf('Bat dau qua trinh lap Newton (Vector & Ma tran):\n');
fprintf('Buoc 0: x = %.10f, y = %.10f\n', X0(1), X0(2));

for i = 0:(N-1)
    % Hiện tại X_old là (x_i, y_i)
    Ji = J(X_old);
    Fi = F(X_old);

    %disp(inv(Ji))

    % 3.1: Kiểm tra điều kiện định thức det(J) != 0
    if det(Ji) == 0
        error('Thong bao: Ma tran Jacobi suy bien (det = 0). Ket thuc thuat toan.');
    end

    % 3.2: Tính X_{i+1} = X_i - J^-1 * F(X_i)
    % Trong Octave, dùng Ji \ Fi thay cho inv(Ji) * Fi để chính xác và nhanh hơn
    X_new = X_old - inv(Ji)*Fi;

    % Cập nhật cho bước sau
    X_old = X_new;

    % In kết quả mỗi bước lặp
    fprintf('Lan lan %d: x = %.10f, y = %.10f\n', i + 1, X_new(1), X_new(2));
end

% B4: In ra giá trị nghiệm gần đúng tìm được
fprintf('\n--- KET QUA CUOI CUNG ---\n');
fprintf('Vec-to nghiem gan dung X = [%.10f; %.10f]\n', X_old(1), X_old(2));
fprintf('Gia tri He ham tai nghiem: f=%.2e, g=%.2e\n', f(X_old(1), X_old(2)), g(X_old(1), X_old(2)));
