clc; clear;
% 1. Định nghĩa ma trận vuông A
A = [4, 1, -1;
     2, 5, -2;
     1, 1,  2];

% 2. Tìm đa thức đặc trưng bằng hàm poly()
% Hàm này trả về các hệ số của đa thức đặc trưng theo thứ tự bậc giảm dần
he_so_da_thuc = poly(A);

% 3. Tìm trị riêng và vectơ riêng bằng hàm eig()
% V: Ma trận chứa các vectơ riêng (mỗi cột là một vectơ riêng)
% D: Ma trận đường chéo chứa các trị riêng tương ứng trên đường chéo chính
[V, D] = eig(A);

% Trích xuất danh sách các trị riêng dạng vectơ từ ma trận đường chéo D
tri_rieng = diag(D);

% 4. Hiển thị kết quả
disp('Ma trận A:');
disp(A);

disp('Các hệ số của đa thức đặc trưng (từ bậc cao đến bậc thấp):');
disp(he_so_da_thuc);

disp('Các trị riêng của ma trận A:');
disp(tri_rieng);

disp('Ma trận vectơ riêng V (mỗi cột là một vectơ riêng tương ứng):');
disp(V);
