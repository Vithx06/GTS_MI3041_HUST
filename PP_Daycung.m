% CHƯƠNG TRÌNH PHƯƠNG PHÁP DÂY CUNG

clc; clear;

% B1: Cấu hình thông số đầu vào
f = @(x) x^5-4096*x-32768;            % Hàm số f(x)
a = 8;                         % Cận trái của khoảng cách ly nghiệm
b = 16;                         % Cận phải của khoảng cách ly nghiệm
eps = 4.2667*10^(-7);               % Sai số cho phép
N0 = 100;                        % Số lần lặp tối đa

% B2: Kiểm tra điều kiện f' và f'' không đổi dấu (Thực hiện ngoài giấy)
% Giả sử điều kiện đã thỏa mãn để tiếp tục thuật toán.

% B3: Khởi tạo x0 và điểm d (điểm cố định)
% Quy tắc: Nếu f(a)*f''(a) < 0 thì x0 = a, d = b. Ngược lại x0 = b, d = a.
x0 = 8;
d = 16;

% B4: Khởi tạo các thông số sai số
i = 0;
x_old = x0;
x = x0;

% B5: Vòng lặp tìm nghiệm
fprintf('\nBat dau qua trinh lap phuong phap day cung:\n');
fprintf(' %-5s   %-20s   %-20s\n', 'k', 'x_k', '|x_k - x_{k-1}|');
fprintf('------------------------------------------------------------\n');

% In dòng k = 0 (Giá trị ban đầu)
fprintf(' %-5d   %-20.10f   %-20s\n', 0, x0, '---');

while i < N0
    % 5.1: Tính x_{i+1} theo công thức dây cung
    tu_so = (d - x_old);
    mau_so = f(d) - f(x_old);
    x_new = x_old - (tu_so / mau_so) * f(x_old);

    % Tính sai số |x_i - x_{i-1}|
    sai_so = abs(x_new - x_old);

    i = i + 1;
    x = x_new; % Cập nhật nghiệm liên tục

    % In dòng dữ liệu của bảng tại vòng lặp i
    fprintf(' %-5d   %-20.10f   %-20.10f\n', i, x_new, sai_so);

    % 5.2: Kiểm tra điều kiện dừng
    if sai_so <= eps
        fprintf('------------------------------------------------------------\n');
        fprintf('=> Thoa man dieu kien sai so tai buoc %d.\n', i);
        break;
    end

    % 5.3: Cập nhật giá trị cũ cho lần lặp kế tiếp
    x_old = x_new;
end

% Kẻ đường kết thúc nếu chạy hết số lần lặp tối đa
if i == N0 && sai_so > eps
    fprintf('------------------------------------------------------------\n');
    fprintf('=> Dung vong lap vi da dat so lan lap toi da N0 = %d.\n', N0);
end

% B6: In ra kết quả cuối cùng
fprintf('\n--- KET QUA CUOI CUNG ---\n');
fprintf('Nghiem gan dung x = %.10f\n', x);
fprintf('So lan lap da thuc hien: %d\n', i);
fprintf('Gia tri f(x) tai nghiem: %.10e\n', f(x));
