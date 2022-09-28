clear all; 
clc; 
A = [2.31 31.49 1.52; 
4.21 22.42 3.85; 
3.49 4.85 28.72]; 
B = [40.95 30.24 42.81]; 
lenA = length(A(:,1)); 
lenB = length(B); 
diag_preob = false; 
for i = 1:lenA 
    for j = 1:lenA 
        sum = 0; 
        if i ~= j 
            sum = sum + abs(A(j,i)); 
        end 
    end 
    if abs(A(i,i)) <= sum 
        diag_preob = true; 
    end 
end 
if ~diag_preob 
x_1 = zeros(1, lenA); 
x_2 = zeros(1, lenA); 
while (max(abs(x_2 - x_1)) >= 0.001 || x_2 - x_1 == 0) 
    x_1 = x_2; 
    for i = 1:lenA 
        x_2(i) = B(i)/A(i,i); 
        for j = 1:lenA 
            if i ~= j 
            x_2(i) = x_2(i) - x_1(j)*A(i,j)/A(i,i); 	  
            end 
        end 
    end 
end 
disp ('Method Yakob:'); 
disp (x_2); 
x_1 = zeros(1, lenA); 
x_2 = zeros(1, lenA); 
while (max(abs(x_2 - x_1)) >= 0.001 || x_2 - x_1 == 0) 
    x_1 = x_2; 
    for i = 1:lenA 
        x_2(i) = B(i)/A(i,i); 
        for j = 1:lenA 
            if i ~= j 
                x_2(i) = x_2(i) - x_2(j)*A(i,j)/A(i,i); 
            end 
        end 
    end 
end 
disp ('Method Gauss-Zeidel:'); 
disp (x_2); 
else 
error ('Matrix does not have the property of diagonal predominance'); 
end 