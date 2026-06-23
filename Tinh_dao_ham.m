% CHƯƠNG TRÌNH TÍNH ĐẠO HÀM (HIỂN THỊ DẠNG PHÂN SỐ TRÊN DƯỚI)
clc; clear;

pkg load symbolic;

% B1: Cấu hình hàm số
syms x;
f = log10(x);

% B2: Tính và rút gọn đạo hàm
df_rut_gon = diff(f, x);

% B3: In tiêu đề
fprintf('\n--- KET QUA DAO HAM f''(x) ---\n');

% B4: Hiển thị dạng phân số "trên dưới" bằng hàm pretty
pretty(df_rut_gon);
