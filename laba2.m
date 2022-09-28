clear all; 
clc; 
Apprxmtn_1 = zeros(2,1); 
Apprxmtn_2 = zeros(2,1); 
Apprxmtn_1_next = zeros(2,1); 
Apprxmtn_2_next = zeros(2,1); 
disp('Initial approximation 1:'); 
for i = 1:2  
    Apprxmtn_1(i) = input(''); 
end 
disp('Initial approximation 2:'); 
for i = 1:2  
    Apprxmtn_2(i) = input(''); 
end 
eps_1 = 0.001; 
eps_2 = 0.001; 
NIT = input('NIT: '); 
disp('No      d1       d2'); 
d1 = 10; 
d2 = 10; 
i = 0; 
while d1 >= eps_1 && d2 >= eps_2 
    i = i + 1; 
    f1 = cos(0.4*Apprxmtn_1(2) + (Apprxmtn_1(1)^2)) + Apprxmtn_1(2)^2 + Apprxmtn_1(1)^2 - 1.6; 
    f2 = 1.5*Apprxmtn_1(1)^2 - (Apprxmtn_1(2)^2)/0.36 - 1; 26         W_1_1 = -
    2*Apprxmtn_1(1)*sin(0.4*Apprxmtn_1(2)+(Apprxmtn_1(1)^2))+2*Apprxmtn_1(1); 
    W_2_1 = 3*Apprxmtn_1(1); 
    W_1_2 = -0.4*sin(0.4*Apprxmtn_1(2)+(Apprxmtn_1(1)^2))+2*Apprxmtn_1(2); 
    W_2_2 = -2*Apprxmtn_1(2)/0.36; 
    discrepancy_vector = [f1; f2]; 
    W = [W_1_1 W_1_2; 
    W_2_1 W_2_2]; 
    W = inv(W); 
    Apprxmtn_1_next = Apprxmtn_1 - W * discrepancy_vector; 
    d1 = max(abs(f1), abs(f2)); 
    d2_v = zeros(2); 
    for j = 1:2 
        if Apprxmtn_1_next(j) < 1 
            d2_v(j) = (Apprxmtn_1_next(j) - Apprxmtn_1(j)); 
        else 
        d2_v(j) = (Apprxmtn_1_next(j) - 
        Apprxmtn_1(j))/Apprxmtn_1_next(j); 
        end 
    end 
    d2 = max(abs(d2_v(1)), abs(d2_v(2))); 
    disp(i+"       "+d1+"      "+d2); 
    Apprxmtn_1 = Apprxmtn_1_next; 
    if i > NIT 
        disp('IER = 2'); 
        break; 
    end 
end 
disp(Apprxmtn_1(1)+" "+Apprxmtn_1(2)); 
d1 = 10; 
d2 = 10; 
i = 0; 
while d1 >= eps_1 && d2 >= eps_2 
    i = i + 1; 
    f1 = cos(0.4*Apprxmtn_2(2) + (Apprxmtn_2(1)^2)) + Apprxmtn_2(2)^2 + Apprxmtn_2(1)^2 - 1.6; 
    f2 = 1.5*Apprxmtn_2(1)^2 - (Apprxmtn_2(2)^2)/0.36 - 1; 60         W_1_1 = -
    2*Apprxmtn_2(1)*sin(0.4*Apprxmtn_2(2)+(Apprxmtn_2(1)^2))+2*Apprxmtn_2(1); 
    W_2_1 = 3*Apprxmtn_2(1); 
    W_1_2 = -0.4*sin(0.4*Apprxmtn_2(2)+(Apprxmtn_2(1)^2))+2*Apprxmtn_2(2); 
    W_2_2 = -2*Apprxmtn_2(2)/0.36; 
    discrepancy_vector = [f1; f2]; 
    W = [W_1_1 W_1_2; 
    W_2_1 W_2_2]; 
    W = inv(W); 
    Apprxmtn_2_next = Apprxmtn_2 - W * discrepancy_vector; 
    d1 = max(abs(f1), abs(f2)); 
    d2_v = zeros(2); 
    for j = 1:2 
        if Apprxmtn_2_next(j) < 1 
            d2_v(j) = (Apprxmtn_2_next(j) - Apprxmtn_2(j)); 
        else 
            d2_v(j) = (Apprxmtn_2_next(j) - 
            Apprxmtn_2(j))/Apprxmtn_2_next(j); 
        end 
    end 
    d2 = max(abs(d2_v(1)), abs(d2_v(2))); 
    disp(i+"       "+d1+"      "+d2); 
    Apprxmtn_2 = Apprxmtn_2_next; 
    if i > NIT 
        disp('IER = 2'); 
        break; 
    end 
end 
disp(Apprxmtn_2(1)+" "+Apprxmtn_2(2)); 