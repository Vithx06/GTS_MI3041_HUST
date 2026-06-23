% CHƯƠNG TRÌNH PHƯƠNG PHÁP LẶP ĐƠN GIẢI HỆ 2 PHƯƠNG TRÌNH
clc; clear;

% B1: Cấu hình thông số đầu vào
phi1 = @(x, y) (sin(x*y))/3; % Hàm số phi1(x, y)
phi2 = @(x, y) (cos(x^2+y^2))/4;       % Hàm số phi2(x, y)

a = 1; b = 2;                            % Khoảng cách ly nghiệm x [a, b]
c = 1; d = 2;                            % Khoảng cách ly nghiệm y [c, d]
x0 = 1/6; y0 = 1/8;                      % Xấp xỉ đầu (x0, y0)
epsilon = 2.482*10^(-2);                  % Sai số cho phép
N0 = 20;                                 % Số lần lặp tối đa

% B2: Hệ số co q
q = 0.2917;

% B3: Khởi tạo
i = 0;
x_old = x0;
y_old = y0;
x = x0; y = y0;

% B4: Vòng lặp tìm nghiệm và in bảng
fprintf('\nBat dau qua trinh lap he phuong trinh:\n');
fprintf(' %-5s   %-15s   %-15s   %-15s\n', 'k', 'x_k', 'y_k', 'Sai so tuong doi');
fprintf('--------------------------------------------------------------------------\n');

% In dòng giá trị ban đầu k = 0
fprintf(' %-5d   %-15.10f   %-15.10f   %-15s\n', 0, x0, y0, '---');

while i < N0
    % 4.1: Tính x_{i+1} và y_{i+1}
    x_new = phi1(x_old, y_old);
    y_new = phi2(x_old, y_old);

    % Tính toán sai số trước khi cập nhật
    tu_so = max(abs(x_new - x_old), abs(y_new - y_old));
    mau_so = max([abs(x_new), abs(y_new)]);
    sai_so_tuong_doi = tu_so / mau_so;

    % Cập nhật nghiệm hiện tại
    x = x_new;
    y = y_new;
    i = i + 1;

    % In dòng dữ liệu của bảng
    fprintf(' %-5d   %-15.10f   %-15.10f   %-15.10f\n', i, x, y, sai_so_tuong_doi);

    % 4.2: Kiểm tra điều kiện dừng
    if sai_so_tuong_doi <= epsilon
        fprintf('--------------------------------------------------------------------------\n');
        fprintf('=> Thoa man dieu kien sai so tai buoc %d.\n', i);
        break;
    end

    % 4.3: Cập nhật giá trị cũ cho lần lặp kế tiếp
    x_old = x_new;
    y_old = y_new;
end

% Kẻ đường kết thúc nếu chạy hết số lần lặp tối đa
if i == N0 && sai_so_tuong_doi > epsilon
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('=> Dung vong lap vi da dat so lan lap toi da N0 = %d.\n', N0);
end

% B5: In ra kết quả cuối cùng
fprintf('\n--- KET QUA CUOI CUNG ---\n');
fprintf('Nghiem gan dung: x = %.10f, y = %.10f\n', x, y);
fprintf('So lan lap da thuc hien: %d\n', i);
