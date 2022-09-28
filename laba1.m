clear all;
clc; 
e = 0.001; 
lyambda = 0; 
lyambda_next = 1000; 
len = input('Matrix size:'); 
disp('Matrix:'); 
for i = 1:len 
    for j = 1:len 
        A(i,j) = input(''); 
    end 
end 
disp('U:'); 
for i = 1:len 
    u_0(i) = input(''); 
end 
u = u_0'; 
u_next = zeros(len); 
while abs(lyambda_next - lyambda) > e 
    lyambda = lyambda_next; 
    u_next = A*u; 
    numerator = dot(u_next,u); 
    denominator = dot(u,u); 
    lyambda_next = numerator/denominator;                                 
    u = u_next; 
    u_next = zeros(len); 
end 
disp('Lyambda_max ='); 
disp(lyambda_next); 
A = inv(A); 
lyambda = 1000;  
lyambda_next = 0; 
u = u_0'; 
u_next = zeros(len); 
while abs(lyambda_next - lyambda) > e 
    lyambda = lyambda_next; 
    u_next = A*u; 
    numerator = dot(u_next,u); 
    denominator = dot(u,u); 
    lyambda_next = numerator/denominator; 
    u = u_next; 
    u_next = zeros(len); 
end 
disp('Lyambda_min =');
disp(lyambda_next); 
