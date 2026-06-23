clc;
clear;

% Cấu hình tham số đầu vào
A = [ 4 2 2;
      2 5 1;
      2 1 6];

n = size(A, 1);
E = 0.00001;          % Sai số hội tụ (epsilon)
Y = [1; 0; 0];   % Vector lặp ban đầu
so_chu_so = 5;      % Số chữ số sau dấu phẩy

% Các hàm phụ

% Chuẩn hóa vector theo chuẩn vô cùng
function V_chuanhoa = Chuan_Hoa_Inf(M)
    [~, idx] = max(abs(M));
    V_chuanhoa = M / M(idx);
end

% Thuật toán xuống thang ma trận để tìm trị riêng tiếp theo
function A_ha_bac = Xuong_Thang(A, VTR, n)
    V = VTR(:);
    theta = eye(n);

    [~, index] = max(abs(V));

    E_col = theta(:, index) - V;
    E_col(index) = 0;
    theta(:, index) = E_col;

    A_ha_bac = theta * A;
end

% Hiển thị lịch sử 3 bước lặp gần nhất trước khi hội tụ
function In_Lich_Su_Lap(B, FMT)
    fprintf('  [LICH SU LAP] 3 Vector cuoi cung truoc khi hoi tu:\n');
    fprintf('  --------------------------------------------------------\n');
    Ma_Tran_Lap = [B{end-2}, B{end-1}, B{end}];
    In_Chuoi = num2str(Ma_Tran_Lap, FMT.fmt_A);

    fprintf('Vector m-2      Vector m-1        Vector m\n');
    disp(In_Chuoi);
    fprintf('  --------------------------------------------------------\n\n');
end

% Tự động tạo chuỗi định dạng hiển thị dựa trên cấu hình số chữ số
FMT.fmt_A = sprintf('     %%.%df     ', so_chu_so);
FMT.fmt_V = sprintf('         %%.%df', so_chu_so);
FMT.fmt_L = sprintf('%%.%df + %%.%dfi', so_chu_so, so_chu_so);
FMT.fmt_P = sprintf('%%%%.%df', so_chu_so);

fprintf('                    MA TRAN BAN DAU A:\n');
fprintf('======================================================================\n');
disp(num2str(A, FMT.fmt_A));
fprintf('======================================================================\n\n');

% Vòng lặp tìm n trị riêng và vector riêng
for i = 1:n
    B = cell(1,1);
    B{1} = Y;
    Z = zeros(n,1);

    % Khởi tạo
    B{2} = Chuan_Hoa_Inf(A * B{1});
    B{3} = Chuan_Hoa_Inf(A * B{2});
    B{4} = Chuan_Hoa_Inf(A * B{3});

    m = 4;
    TH = 4; % Không thỏa mãn trường hợp nào thì là TH4

    while true

        sai_so_tron = 10^so_chu_so;
        F = round(B{m} * sai_so_tron) / sai_so_tron;

        if all(F == Z)
            error('Vector hoi tu ve 0, thuat toan dung.');
        end

        % TH1: B(m) ≈ B(m-1)
        A_diff1 = B{m} - B{m-1};
        if max(abs(A_diff1)) <= E
            TH = 1;
            break;
        end

        % TH3: B(m) ≈ B(m-2)
        A_diff2 = B{m} - B{m-2};
        if max(abs(A_diff2)) <= E
            TH = 3;
            break;
        end

        B{m+1} = Chuan_Hoa_Inf(A * B{m});

        % Cập nhật vector mới nhất
        m = m + 1;

    end

    m_iters = m;

% Xử lý từng trường hợp
    if TH == 1
        fprintf('======================================================================\n');
        fprintf('=> Giai lamda_%d thuoc TH %d (Tri rieng thuc troi)\n', i, TH);
        fprintf('======================================================================\n');

        In_Lich_Su_Lap(B, FMT);

        VTR = B{end};
        AX = A * VTR;

        % Tính trị riêng bằng thương số tại vị trí có giá trị lớn nhất
        [~, idx] = max(abs(VTR));
        lamda = AX(idx) / VTR(idx);

        fprintf('  * Ket qua tai buoc lap thu %d:\n', m_iters - 1);
        fprintf(['    + lamda_%d = %.' num2str(so_chu_so) 'f\n'], i, lamda);
        fprintf('    + Vector rieng tuong ung:\n');
        disp(num2str(VTR, FMT.fmt_V));

        % Xuống thang ma trận để loại bỏ trị riêng vừa tìm được
        A = Xuong_Thang(A, VTR, n);

        fprintf('  * Ma tran sau khi Xuong thang:\n');
        disp(num2str(A, FMT.fmt_A));
        fprintf('______________________________________________________________________\n\n');

        % Đổi vector ban đầu ngẫu nhiên để tránh suy biến ở vòng lặp sau
        Y = rand(n, 1);

    elseif TH == 3
        fprintf('======================================================================\n');
        fprintf('=> Giai lamda_%d thuoc TH %d (Hai tri rieng doi nhau)\n', i, TH);
        fprintf('======================================================================\n');

        In_Lich_Su_Lap(B, FMT);

        AmY = A * B{end};
        Am1Y = A * AmY;
        Am2Y = A * Am1Y;

        [~, idx] = max(abs(AmY));

        % Tính trị riêng cho trường hợp đối nhau: lambda = sqrt(A^(m+2)Y / A^mY)
        lamda_N = sqrt(Am2Y(idx) / AmY(idx));

        % Tính vector riêng tương ứng
        VTR_N = Am1Y + lamda_N * AmY;
        VTR_N = Chuan_Hoa_Inf(VTR_N);

        fprintf('  * Ket qua tai buoc lap thu %d:\n', m_iters - 1);
        fprintf(['    + lamda_%d = %.' num2str(so_chu_so) 'f\n'], i, lamda_N);
        fprintf('    + Vector rieng tuong ung:\n');
        disp(num2str(VTR_N, FMT.fmt_V));

        % Xuống thang ma trận
        A = Xuong_Thang(A, VTR_N, n);

        fprintf('  * Ma tran sau khi Xuong thang:\n');
        disp(num2str(A, FMT.fmt_A));
        fprintf('______________________________________________________________________\n\n');

        Y = rand(n, 1);

    elseif TH == 4
        fprintf('======================================================================\n');
        fprintf('=> Giai lamda_%d thuoc TH %d (Cap tri rieng phuc lien hop)\n', i, TH);
        fprintf('======================================================================\n\n');

        M = B{end};
        AmY = A * M;
        Am1Y = A * AmY;

        % Lấy các thành phần của vector để thiết lập hệ phương trình sai phân
        a_1 = Am1Y(1); a_2 = Am1Y(2);
        b_1 = AmY(1);  b_2 = AmY(2);
        c_1 = M(1);    c_2 = M(2);

        % Tính các hệ số của phương trình đặc trưng bậc 2: Z^2 + b_coef*Z + c_coef = 0
        a_coef = 1;
        b_coef = (a_1 * c_2 - c_1 * a_2) / (c_1 * b_2 - b_1 * c_2);
        c_coef = (b_1 * a_2 - a_1 * b_2) / (c_1 * b_2 - b_1 * c_2);

        fmt_ptdt = sprintf('  * Phuong trinh dac trung dactinh: Z^2 + (%s)Z + (%s) = 0 \n\n', FMT.fmt_P, FMT.fmt_P);
        fprintf(fmt_ptdt, b_coef, c_coef);

        denta = b_coef^2 - 4 * a_coef * c_coef;

        if denta >= 0
            disp('  * Phuong trinh co denta >= 0 (Khong phai nghiem phuc). Dung thuat toan.');
        else
            % Tìm 2 trị riêng phức liên hợp
            lamda_1 = complex(-b_coef / (2 * a_coef), -sqrt(abs(denta)) / (2 * a_coef));
            lamda_2 = complex(-b_coef / (2 * a_coef), sqrt(abs(denta)) / (2 * a_coef));

            % Tìm 2 vector riêng phức tương ứng
            VTR_1 = Am1Y - lamda_1 * AmY;
            VTR_2 = Am1Y - lamda_2 * AmY;

            fprintf(['    + lamda_%d = ' FMT.fmt_L '\n'], i, real(lamda_1), imag(lamda_1));
            fprintf('    + Vector rieng tuong ung:\n');
            disp(num2str(VTR_1, FMT.fmt_V));

            fprintf('    %s\n', repmat('-', 1, 50));

            fprintf(['    + lamda_%d = ' FMT.fmt_L '\n'], i + 1, real(lamda_2), imag(lamda_2));
            fprintf('    + Vector rieng tuong ung:\n');
            disp(num2str(VTR_2, FMT.fmt_V));
            fprintf('======================================================================\n\n');
        end
        break; % Kết thúc vòng lặp vì TH4 giải quyết cùng lúc 2 trị riêng cuối cùng
    end
end
