% 1. Định nghĩa ma trận hệ số A và vectơ kết quả B
%

clc;
clear;

A = [
    33 -1 2 1 -3 2;
    -5 33 -1 5 2 -3;
    1 3 33 2 2 1;
    2 1 -1 33 0 4;
    3 -8 4 3 33 2;
    2 3 1 4 2 33
];

B = [
    29.606 24.8297;
    3.5764 8.5379;
    26.9654 11.342;
    34.3139 31.1558;
    16.2633 15.5075;
    19.7982 29.9398
];

% 2. Giải hệ phương trình bằng toán tử chia trái (Khuyên dùng)
X1 = A \ B;

% 3. Giải hệ phương trình bằng hàm linsolve có sẵn
X2 = linsolve(A, B);

% 4. Hiển thị kết quả
disp('Kết quả giải bằng toán tử A \ B:');
disp(X1);

disp('Kết quả giải bằng hàm linsolve(A, B):');
disp(X2);
