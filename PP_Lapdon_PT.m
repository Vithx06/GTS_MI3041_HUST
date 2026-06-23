% CHƯƠNG TRÌNH PHƯƠNG PHÁP LẶP ĐƠN PT

clc; clear;

% B1: Cấu hình thông số đầu vào
phi = @(x) 0.1*exp(sqrt(-x))-3;        % Hàm số phi(x)
a = -3;                                % Cận trái
b = -2;                                % Cận phải
x0 = -2;                               % Giá trị xấp xỉ đầu
epsilon = 10^(-10);                     % Sai số cho phép
N0 = 100;                                % Số lần lặp tối đa

% B2: Tính hệ số co q: Tính tay
q = 0.163;

% B3: Kiểm tra điều kiện hội tụ
if q >= 1
    fprintf('Thong bao: Khong thoa man DK phuong phap lap don (q >= 1).\n');
    return;
end

% B4: Khởi tạo
i = 0;
x_old = x0;
eps_hat = epsilon * (1 - q) / q;

% B5: Vòng lặp tìm nghiệm
fprintf('\nBat dau qua trinh lap:\n');
fprintf(' %-5s   %-20s   %-20s\n', 'k', 'x_k', '|x_k - x_{k-1}|');
fprintf('-----------------------------------------------------------\n');

% In dòng k = 0
fprintf(' %-5d   %-20.10f   %-20s\n', 0, x0, '---');

while i < N0
    % 5.1: Tính x_{i+1}
    x_new = phi(x_old);
    i = i + 1;

    % Tính sai số |x_i - x_{i-1}|
    sai_so = abs(x_new - x_old);

    % In dòng dữ liệu (Căn lề trái để dễ nhìn)
    fprintf(' %-5d   %-20.10f   %-20.10f\n', i, x_new, sai_so);

    % Cập nhật biến x để lưu kết quả cuối
    x = x_new;

    % 5.2: Kiểm tra điều kiện dừng
    if sai_so <= eps_hat
        fprintf('-----------------------------------------------------------\n');
        fprintf('=> Thoa man dieu kien sai so tai buoc %d.\n', i);
        break;
    end

    % 5.4: Cập nhật giá trị cũ cho lần lặp kế tiếp
    x_old = x_new;
end

% Kẻ đường kết thúc nếu đạt tối đa lần lặp
if i == N0 && sai_so > eps_hat
    fprintf('-----------------------------------------------------------\n');
    fprintf('=> Dung vong lap vi da dat so lan lap toi da N0 = %d.\n', N0);
end

% B6: In ra kết quả cuối cùng
fprintf('\n--- KET QUA CUOI CUNG ---\n');
fprintf('Nghiem gan dung x = %.10f\n', x);
fprintf('So lan lap da thuc hien: %d\n', i);
fprintf('Gia tri phi(x) tai nghiem: %.10f\n', phi(x));
