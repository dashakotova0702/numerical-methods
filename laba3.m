clear all; 
clc; 
A = [2.31 31.49 1.52; 
4.21 22.42 3.85; 
3.49 4.85 28.72]; 
B = [40.95; 30.24; 42.81]; 
A_norm = A; 
B_norm = B; 
x = 0; 
if (det(A) == 0) 
error('det A == 0');  
end 
lenA = length(A(:,1)); 
lenB = length(B); 
if lenB ~= lenA 
error('length B != length A'); 
end 
for i = 1:lenA 
    for j = 1:lenA 
        if j - i <= 0 
            continue;  	 	 	 	 	 	 	 	 	 
        end 
        a = A(i,i); 
        b = A(j,i); 
        for k = 1:lenA 
            A(j,k) = a*(A(j,k))-b*(A(i,k)); 
        end  
        B(j,1) = a*(B(j,1))-b*(B(i,1)); 
    end 
end 
for i = lenA:-1:1 
    x(i) = B(i, 1); 
    for j = [1:i-1, i+1:lenA] 
        x(i) = x(i) - A(i, j) * x(j); 
    end 
    x(i) = x(i)/A(i,i); 
end 
F = zeros(1,3); 
for i = 1:lenA 
    for j = 1:lenA 
        F(i) = F(i) + A_norm(i,j) * x(j); 
    end 
    F(i) = F(i) - B_norm(i,1); 
end 
norm = max(F); 
disp('x:');
disp(x); 
disp('F:'); 
disp(F); 
disp('Norm:'); 
disp(norm); 