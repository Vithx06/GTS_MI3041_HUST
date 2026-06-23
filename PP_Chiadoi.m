% CHƯƠNG TRÌNH PHƯƠNG PHÁP CHIA ĐÔI - XUẤT BẢNG
clc; clear;

% ================= B1: CẤU HÌNH THÔNG SỐ =================
f = @(x) 2^x - cos(x + 2);  % Hàm số f(x)
a = -2;                     % Cận trái ban đầu
b = -0.5;                   % Cận phải ban đầu
epsilon = 1e-10;            % Sai số (để dừng nếu cần)
N0 = 8;                     % Số lần lặp tối đa (N)

% ================= B2: KIỂM TRA ĐIỀU KIỆN =================
fa = f(a);
fb = f(b);

if fa * fb > 0
    fprintf('Thông báo: f(a) và f(b) cùng dấu. Không thể thực hiện.\n');
    return;
end

% ================= B3: IN TIÊU ĐỀ BẢNG =================
fprintf(' %-3s | %-12s | %-12s | %-12s | %-10s\n', 'n', 'an', 'bn', 'cn', 'sign(f(cn))');
fprintf('----------------------------------------------------------------------\n');

% Khởi tạo biến lưu giá trị cũ để so sánh
a_old = NaN;
b_old = NaN;

% ================= B4: VÒNG LẶP VÀ IN KẾT QUẢ =================
for n = 0:N0
    % Tính điểm giữa cn
    cn = (a + b) / 2;
    fcn = f(cn);

    % Xử lý hiển thị cho an
    if n > 0 && a == a_old
        str_a = '-----';
    else
        str_a = sprintf('%.8f', a);
    end

    % Xử lý hiển thị cho bn
    if n > 0 && b == b_old
        str_b = '-----';
    else
        str_b = sprintf('%.8f', b);
    end

    % In dòng kết quả hiện tại
    fprintf(' %-3d | %-12s | %-12s | %-12.8f | %-10d\n', n, str_a, str_b, cn, sign(fcn));

    % Lưu giá trị hiện tại trước khi cập nhật cho bước sau
    a_old = a;
    b_old = b;

    % Cập nhật khoảng cách ly nghiệm cho bước n + 1
    if f(a) * fcn < 0
        b = cn;
    else
        a = cn;
    end

    % Kiểm tra điều kiện dừng sớm (nếu cần)
    if abs(fcn) < eps || (b - a) / 2 < epsilon
        % Nếu muốn dừng khi đạt sai số thì thêm break ở đây
    end
end

fprintf('----------------------------------------------------------------------\n');
fprintf('Nghiem gan dung cuoi cung: x = %.10f\n', cn);
