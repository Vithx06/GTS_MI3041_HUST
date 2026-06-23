% CHƯƠNG TRÌNH PHƯƠNG PHÁP TIẾP TUYẾN (NEWTON)
clc; clear;

% B1: Cấu hình thông số đầu vào
f = @(x) 1.2*x^5-2.57*x+2;           % Hàm số f(x)
df = @(x) 6*x^4-2.57;                % BỔ SUNG: Đạo hàm f'(x) để tính toán
a = -1.5;                         % Cận trái
b = -1;                         % Cận phải
eps = 9.20279*10^(-5);             % Sai số cho phép
N0 = 100;                       % Số lần lặp tối đa

% B3: Kiểm tra điều kiện Fourier để chọn x0
x0 = -1.5;

% B4: Khởi tạo các thông số sai số
i = 0;
x_old = x0;
x = x0;

% B5: Vòng lặp tìm nghiệm
fprintf('\nBat dau qua trinh lap phuong phap tiep tuyen:\n');
fprintf(' %-5s   %-20s   %-20s\n', 'k', 'x_k', '|x_k - x_{k-1}|');
fprintf('------------------------------------------------------------\n');

% In dòng k = 0 (Giá trị khởi đầu)
fprintf(' %-5d   %-20.10f   %-20s\n', 0, x0, '---');

while i < N0
    % 5.1: Tính x_{i+1} theo công thức Newton
    x_new = x_old - f(x_old) / df(x_old);

    % Tính sai số |x_i - x_{i-1}|
    sai_so = abs(x_new - x_old);

    i = i + 1;
    x = x_new; % CẬP NHẬT NGHIỆM LIÊN TỤC

    % In dòng dữ liệu của bảng tại vòng lặp i
    fprintf(' %-5d   %-20.10f   %-20.10f\n', i, x_new, sai_so);

    % 5.2: Kiểm tra điều kiện dừng
    if sai_so <= eps
        fprintf('------------------------------------------------------------\n');
        fprintf('=> Thoa man dieu kien sai so tai buoc %d.\n', i);
        break;
    end

    % 5.3: Cập nhật giá trị cho lần lặp kế tiếp
    x_old = x_new;
end

% Kẻ đường kết thúc nếu chạy hết số lần lặp tối đa mà chưa dừng
if i == N0 && sai_so > eps
    fprintf('------------------------------------------------------------\n');
    fprintf('=> Dung vong lap vi da dat so lan lap toi da N0 = %d.\n', N0);
end

% B6: In ra kết quả cuối cùng
fprintf('\n--- KET QUA CUOI CUNG ---\n');
fprintf('Nghiem gan dung x = %.10f\n', x);
fprintf('So lan lap da thuc hien: %d\n', i);
fprintf('Gia tri f(x) tai nghiem: %.10e\n', f(x));
