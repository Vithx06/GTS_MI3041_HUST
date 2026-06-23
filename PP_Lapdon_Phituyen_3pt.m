% CHƯƠNG TRÌNH PHƯƠNG PHÁP LẶP ĐƠN GIẢI HỆ 3 PHƯƠNG TRÌNH
clc; clear;

% B1: Cấu hình thông số đầu vào
phi1 = @(x, y, z) (cos(y + z))/5;
phi2 = @(x, y, z) (sin(x + z))/6;
phi3 = @(x, y, z) (exp(-(x + y)))/8;

a = 0; b = 1;
c = 0; d = 1;
e = 0; f = 1;

x0 = 0.1; y0 = 0.1; z0 = 0.1;
epsilon = 10^(-5);
N0 = 20;

% B2: Hệ số co q
q = 0.35;

% B3: Khởi tạo
i = 0;
x_old = x0;
y_old = y0;
z_old = z0;

x = x0; y = y0; z = z0;

% B4: Vòng lặp tìm nghiệm và in bảng
fprintf('\nBat dau qua trinh lap he 3 phuong trinh:\n');
fprintf(' %-5s   %-15s   %-15s   %-15s   %-15s\n', 'k', 'x_k', 'y_k', 'z_k', 'Sai so tuong doi');
fprintf('------------------------------------------------------------------------------------------\n');

% In dòng giá trị ban đầu k = 0
fprintf(' %-5d   %-15.10f   %-15.10f   %-15.10f   %-15s\n', 0, x0, y0, z0, '---');

while i < N0
    % 4.1: Tính x_{i+1}, y_{i+1}, z_{i+1}
    x_new = phi1(x_old, y_old, z_old);
    y_new = phi2(x_old, y_old, z_old);
    z_new = phi3(x_old, y_old, z_old);

    % Tính toán sai số tương đối
    tu_so = max([abs(x_new - x_old), abs(y_new - y_old), abs(z_new - z_old)]);
    mau_so = max([abs(x_new), abs(y_new), abs(z_new)]);
    sai_so_tuong_doi = tu_so / mau_so;

    % Cập nhật nghiệm hiện tại
    x = x_new;
    y = y_new;
    z = z_new;
    i = i + 1;

    % In dòng dữ liệu của bảng
    fprintf(' %-5d   %-15.10f   %-15.10f   %-15.10f   %-15.10f\n', i, x, y, z, sai_so_tuong_doi);

    % 4.2: Kiểm tra điều kiện dừng
    if sai_so_tuong_doi <= epsilon
        fprintf('------------------------------------------------------------------------------------------\n');
        fprintf('=> Thoa man dieu kien sai so tai buoc %d.\n', i);
        break;
    end

    % 4.3: Cập nhật giá trị cũ cho lần lặp kế tiếp
    x_old = x_new;
    y_old = y_new;
    z_old = z_new;
end

% Kẻ đường kết thúc nếu chạy hết số lần lặp tối đa
if i == N0 && sai_so_tuong_doi > epsilon
    fprintf('------------------------------------------------------------------------------------------\n');
    fprintf('=> Dung vong lap vi da dat so lan lap toi da N0 = %d.\n', N0);
end

% B5: In ra kết quả cuối cùng
fprintf('\n--- KET QUA CUOI CUNG ---\n');
fprintf('Nghiem gan dung:\n');
fprintf('x = %.10f\n', x);
fprintf('y = %.10f\n', y);
fprintf('z = %.10f\n', z);
fprintf('So lan lap da thuc hien: %d\n', i);
