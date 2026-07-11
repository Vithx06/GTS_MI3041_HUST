% Xoa man hinh va bo nho
clc; clear all;

% ==========================================
% 1. CAU HINH BAI TOAN (Dau vao)
% ==========================================
% Ban co the thay doi cac gia tri duoi day de giai bai toan cua minh
A = [1.05, 0.05, 0.02;
     0.02, 1.04, 0.01;
     0.01, 0.02, 1.03];

B = [1; 2; 3];

X0 = [0; 0; 0];

p = inf; % Chuan p (co the la 1, 2, hoac inf)
epsilon = 1e-6; % Sai so cho phep
N0 = 100; % So buoc lap toi da

% In cau hinh ban dau
printf("--- CAU HINH DAU VAO ---\n");
printf("Ma tran A:\n"); disp(num2str(A, '%.7f  '));
printf("Ma tran B:\n"); disp(num2str(B, '%.7f  '));
printf("Xap xi dau X(0):\n"); disp(num2str(X0, '%.7f  '));
printf("Chuan p: %s\n", num2str(p));
printf("Sai so epsilon: %.7f\n", epsilon);
printf("So buoc lap toi da N0: %d\n\n", N0);

% ==========================================
% 2. THUAT TOAN LAP DON
% ==========================================
n = size(A, 1);
I = eye(n);

% Buoc 2: Tinh ma tran lap C va he so co q
C = I - A;
q = norm(C, p);

printf("--- BUOC 2: TINH MA TRAN LAP VA HE SO CO ---\n");
printf("Ma tran lap C = I - A:\n"); disp(num2str(C, '%.7f  '));
printf("He so co q = ||C||_p = %.7f\n\n", q);

% Buoc 3: Kiem tra dieu kien co
if q >= 1 || q == 0
    printf("THONG BAO: q = %.7f khong thoa man dieu kien co (0 < q < 1).\n", q);
    printf("Thuat toan ket thuc som.\n");
    return;
end

% Buoc 4: Khoi tao bien dem va xap xi
k = 1;
X_prev = X0;

printf("--- BAT DAU QUA TRINH LAP ---\n");

% Buoc 5: Vong lap chinh
while k <= N0
    % 5.1: Tinh xap xi tiep theo X^(k)
    X_curr = C * X_prev + B;

    % 5.2: Tinh sai so eps tai buoc k
    eps_k = (q / (1 - q)) * norm(X_curr - X_prev, p);

    % In ket qua tung buoc
    printf(">> Buoc k = %d:\n", k);
    printf("X^(%d) =\n", k);
    disp(num2str(X_curr, '%.7f  '));
    printf("Sai so eps = %.7f\n\n", eps_k);

    % 5.3: Kiem tra dieu kien dung
    if eps_k <= epsilon
        printf("=> Thoa man dieu kien dung tai buoc %d.\n\n", k);
        break;
    end

    % 5.4: Chuan bi cho buoc lap tiep theo
    X_prev = X_curr;
    k = k + 1;
end

% Buoc 6: Kiem tra loi khong hoi tu (Vuot qua N0)
if k > N0
    printf("THONG BAO: Thuat toan khong hoi tu sau %d buoc. Ket thuc!\n", N0);
    return;
end

% Buoc 7: Tong hop va in ket qua cuoi cung
printf("==========================================\n");
printf("             TONG HOP KET QUA             \n");
printf("==========================================\n");
printf("- So buoc lap k: %d\n", k);
printf("- Sai so dat duoc eps: %.7f\n", eps_k);
printf("- Nghiem gan dung X^(%d) la:\n", k);
disp(num2str(X_curr, '%.7f  '));
printf("==========================================\n");
