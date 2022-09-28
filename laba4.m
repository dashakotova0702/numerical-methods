clear all;
clc;
disp('Lx:');
Lx = input('');
disp('Ly:');
Ly = input('');
disp('T:');
T = input('');
disp('c:');
c = input('');
ft = "100*sin(10*t)";
f_t = inline(ft, 't');
explicit_method(Lx, Ly, T, c, f_t);
implicit_method(Lx, Ly, T, c, f_t);
function explicit_method(Lx, Ly, T, c, f_t)
    delta_x = 0.1;
    delta_y = 0.1;
    delta_t = (delta_x^2/c)*0.1;
    length_x = Lx/delta_x + 1;
    length_y = Ly/delta_y + 1;
    x = 0:delta_x:Lx;
    y = 0:delta_y:Ly;
    u = zeros(length_x, length_y);
    for i = 1:length_x
        for j = 1:length_y
            u(i, j) = 30;
        end
    end
    [X,Y] = meshgrid(y,x); 
    gr = surf(X, Y, u);
    pause(delta_t);
    delete(gr);
    t = 0;
    while t <= T
        for i = 1:length_x
            for j = 1:length_y
                if i == 1
                    u_x(i, j) = (u(i+1, j)-2*(u(i,j))+30)/delta_x^2;
                else
                    if i == length_x
                        u_x(i, j) = (100-2*(u(i,j))+u(i-1, j))/delta_x^2;
                    else
                        u_x(i, j) = (u(i+1, j)-2*(u(i,j))+u(i-1, j))/delta_x^2;
                    end
                end
                if j == 1
                    u_y(i, j) = (u(i, j+1)-u(i,j))/delta_y^2;
                else
                    if j == length_y
                        u_y(i, j) = (-u(i,j)+u(i, j-1))/delta_y^2;
                    else
                        u_y(i, j) = (u(i, j+1)-2*(u(i,j))+u(i, j-1))/delta_y^2;
                    end
                end
            end
        end
        u = u + delta_t*(c .* (u_x + u_y) + f_t(t));
        gr = surf(X, Y, u);
        pause(delta_t/3);
        delete(gr);
        disp(t);
        t = t + delta_t;
    end
    close;
end
function implicit_method(Lx, Ly, T, c, f_t)
    delta_x = 0.1;
    delta_y = 0.1;
    delta_t = (delta_x^2/c)*0.1;
    length_x = Lx/delta_x + 1;
    length_y = Ly/delta_y + 1;
    x = 0:delta_x:Lx;
    y = 0:delta_y:Ly;
    u = zeros(length_x, length_y);
    for i = 1:length_x
        for j = 1:length_y
            u(i, j) = 30;
        end
    end
    [X,Y] = meshgrid(y,x); 
    gr = surf(X, Y, u);
    pause(delta_t);
    delete(gr);
    t = 0;
    while t <= T
        u_coeff = zeros((length_x),(length_x));
        for i = 1:length_x
            if i == 1
                    u_coeff(i, i) = 1+2*c*delta_t/delta_x^2;
                    u_coeff(i, i+1) = -c*delta_t/delta_x^2;
                    u_x_prev(i) = u(i,1) + delta_t*(c*30/delta_x^2 - f_t(t+delta_t));
            else
                if i == length_x
                    u_coeff(i, i-1) = -c*delta_t/delta_x^2;
                    u_coeff(i, i) = 1+2*c*delta_t/delta_x^2;
                    u_x_prev(i) = u(i,1) + delta_t*(c*100/delta_x^2 - f_t(t+delta_t));
                else
                    u_coeff(i, i-1) = -c*delta_t/delta_x^2;
                    u_coeff(i, i) = 1+2*c*delta_t/delta_x^2;
                    u_coeff(i, i+1) = -c*delta_t/delta_x^2;
                    u_x_prev(i) = u(i,1) - delta_t*f_t(t+delta_t);
                end
            end
        end
        u_x_next = u_coeff\u_x_prev';
        for j = 1:length_y
            u(:,j) = u_x_next;
        end
        gr = surf(X, Y, u);
        pause(delta_t/3);
        delete(gr);
        disp(t);
        t = t + delta_t;
    end
    close;
end
