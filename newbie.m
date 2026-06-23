#Phương pháp chia đôi giải f(x) = 0
%{
Các điều kiện để thực hiện phương pháp ( giả sử đã thỏa mãn ):
1. (a,b) là khoảng cách ly nghiệm
2. f(x) liên tục trên (a,b)
3. f(a).f(b) < 0
%}

% A. Cấu hình các tham số : hàm số - khoảng cách ly (a,b) - epxilon hoặc n - số chữ số sau dấu ,
f = @(x) 4*exp(2*x) - 3*x^2 + 7; a = -2; b = -1;
eps = 0.5*10^(-2) ;       % TH1: Nhập epsilon (Nếu dùng n thì để = 0)
n_fix = 0 ;     % TH2: Nhập số bước n (Nếu dùng eps thì để = 0)
prec = 12; cw = prec + 3; % Độ rộng cột

% --- PHẢI TÍNH fa, fb Ở ĐÂY ĐỂ CÓ DẤU CHO TIÊU ĐỀ ---
fa = f(a); fb = f(b);
if fa * fb >= 0; printf("\n[!] f(a) va f(b) cung dau. Khong the khep nghiem!\n"); return; end

% --- TÍNH n_total (Dùng if-else thuần túy để tránh lỗi eps=0) ---
if n_fix > 0; n_total = n_fix; m_txt = "Co dinh n"; else; n_total = ceil(log2(abs(b-a)/eps)); m_txt = "Theo eps"; end

% ==========================================================
% B. IN TIÊU ĐỀ
% ==========================================================
sa = {"(-)", "(+)"}{(fa > 0) + 1}; sb = {"(-)", "(+)"}{(fb > 0) + 1};
printf("\nChe do: %s (n=%d)\n", m_txt, n_total);
printf("%-3s | %-*s | %-*s | %-*s | %-5s\n", "n", cw, ["an", sa], cw, ["bn", sb], cw, "cn", "f(cn)");
disp(repmat('-', 1, 20 + 3 * cw));

% ==========================================================
% C. THỰC HIỆN THUẬT TOÁN
% ==========================================================
ap = bp = -9e9;

for n = 0 : n_total
    % Mẹo Indexing: {Gia_tri_khi_SAI, Gia_tri_khi_DUNG}{(Dieu_kien) + 1}
    str_a = {sprintf("%.*f", prec, a), repmat('-', 1, cw-2)}{(abs(a-ap) < 1e-15) + 1};
    str_b = {sprintf("%.*f", prec, b), repmat('-', 1, cw-2)}{(abs(b-bp) < 1e-15) + 1};
    c = (a + b) / 2; z = f(c);

    if n == n_total % Dòng cuối cùng: cn và f(cn) để dấu \
        printf("%-3d | %-*s | %-*s | %-*s | %-5s\n", n, cw, str_a, cw, str_b, cw, "\\", "\\");
    else            % Dòng bình thường: xác định dấu f(c) bằng mẹo Indexing
        printf("%-3d | %-*s | %-*s | %-*.*f | %-5s\n", n, cw, str_a, cw, str_b, cw, prec, c, {"-", "+"}{(z > 0) + 1});
    end

    ap = a; bp = b; % Lưu giá trị cũ để dòng sau so sánh
    if z * f(ap) < 0; b = c; else; a = c; end % Cập nhật đầu mút
end

disp(repmat('-', 1, 20 + 3 * cw));
printf("Nghiem xap xi: %.*f\n", prec, b );
