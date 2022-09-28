clear all; 
clc; 
D = 0:0.1:0.9; 
v = [0.957, 0.969, 0.976, 0.978, 0.975, 0.968, 0.954, 0.939, 0.918, 0.894]; 
N = 10; 
m = 2; 
POWERX = zeros(2*m); 
for k = 1:2*m 
    for i = 1:N 
        POWERX(k) =  POWERX(k) + D(i).^k; 
    end 
end 
SUMX = zeros(m+1, m+1); 
for l = 1:m+1 
    for j = 1:m+1 
        if (j == 1 && l == 1) 
            SUMX(l, j) = N; 
        else 
            SUMX(l, j) = POWERX(l + j - 2);  	 	 	 	 	 	
        end 
    end 
end 
PRAW = zeros(m+1); 
for l = 1:m+1  
    for i = 1:N 
        PRAW(l) = PRAW(l) + v(i).*(D(i).^(l-1)); 
    end 
end 
coeff = SUMX.'\PRAW; 
coeff = coeff.'; 
S2 = 1; 
for i = 1:N 
    S2_ = v(i); 
    for j = 1:m+1 
        S2_ = S2_ - coeff(j).*(D.^(j-1));  
    end 
    S2 = S2 + (S2_.^2); 
end 
S2 = S2./(N-m-1); 
sigma = sqrt(S2); 
disp(S2); 
disp(sigma); 
x = 0:0.01:1; 
y = coeff(1, 3)*x.^2 + coeff(1, 2)*x + coeff(1, 1); 
figure; 
plot (D, v, '-o', 'Color', 'r'); 
hold on 
plot (x, y, 'b'); 