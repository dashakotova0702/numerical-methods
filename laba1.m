clear all; 
clc; 
hold on 
x = -4*pi:0.001:4*pi; 
y1 = sin(x); 
y = -4*pi:0.001:4*pi; 
plot(x, y1); 
plot(x, y); 
figure; 
n = 1; 
hold off 
while (abs(y(length(y))-y1(length(y))) >= 0.001)  	 	 	 	
    hold on 
    plot(x, y1); 
    y = y + (((-1).^n).*(x.^(2*n+1)))/(factorial(2*n+1)); 
    plot (x, y); 
    figure; 
    xlabel('x[]') 
    ylabel('y[]') 
    n = n + 1; 
    disp(abs(y(length(y))-y1(length(y)))); 
    hold off 
end 