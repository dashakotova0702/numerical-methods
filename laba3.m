clear all; 	 
clc; 
explicit_method (2, [0, -0.412], 0.1, 1, 0.1, 0.001, ["-u1*u2+sin(t)/t", "-(u2.^2)+3.5*t/(1+t.^2)"]); 
figure; 
implicit_method (2, [0, -0.412], 0.005, 0.001, 1, 0.001, ["-u1*u2+sin(t)/t", "-(u2.^2)+3.5*t/(1+t.^2)"]); 
figure; 
shihman_method (2, [0, -0.412], 0.005, 0.001, 1, 0.001, ["-u1*u2+sin(t)/t", "-(u2.^2)+3.5*t/(1+t.^2)"]); 
figure; 
explicit_method (2, [1, 0], 0.01, 1, 0.01, 0.001, ["u2-(2.25*u1+u2)*u1", "exp(u1)-(u1+2.25*u2)*u1"]); 
figure; 
implicit_method (2, [1, 0], 0.005, 0.001, 1, 0.001, ["u2-(2.25*u1+u2)*u1", "exp(u1)-(u1+2.25*u2)*u1"]); 
figure; 
shihman_method (2, [1, 0], 0.005, 0.001, 1, 0.001, ["u2-(2.25*u1+u2)*u1", "exp(u1)-(u1+2.25*u2)*u1"]); 
figure; 
explicit_method (3, [1, 1, 1], 0.01, 1, 0.01, 0.001, ["(2.25-1.25)/1.25*u2*u3", "(1.25+2.25)/2.25*u1*u3", "(1.25-2.25)/1.25*u1*u2"]); 
figure; 
implicit_method (3, [1, 1, 1], 0.005, 0.001, 1, 0.001, ["(2.25-1.25)/1.25*u2*u3", "(1.25+2.25)/2.25*u1*u3", "(1.25-2.25)/1.25*u1*u2"]); 
figure; 
shihman_method (3, [1, 1, 1], 0.005, 0.001, 1, 0.001, ["(2.25-1.25)/1.25*u2*u3", "(1.25+2.25)/2.25*u1*u3", "(1.25-2.25)/1.25*u1*u2"]); 


function x_next = Newton (count, x0, vector_f, t_now) 
    syms u1 u2 u3 t; 
    c = 0; 
    if count == 3 
        x_next = [1000, 1000, 1000]; 
        while abs(x_next(1)-x0(1)) > 0.01 && abs(x_next(2)-x0(2)) > 0.01 && abs(x_next(3)-x0(3)) > 0.01 
            c = c + 1; 
            if c > 5 
                return; 
            end 
            x_next = x0; 
            x0 = x_next; 
            f0 = [double(subs((subs((subs(vector_f{1}, u1, x0(1))), u2, x0(2))), u3, x0(3)));  
            double(subs((subs((subs(vector_f{2}, u1, x0(1))), u2, x0(2))), u3, x0(3)));  
            double(subs((subs((subs(vector_f{3}, u1, x0(1))), u2, x0(2))), u3, x0(3)))]; 
            df1u1 = diff(vector_f{1}, u1); 
            df1u2 = diff(vector_f{1}, u2); 
            df1u3 = diff(vector_f{1}, u3); 
            df2u1 = diff(vector_f{2}, u1); 
            df2u2 = diff(vector_f{2}, u2); 
            df2u3 = diff(vector_f{2}, u3); 
            df3u1 = diff(vector_f{3}, u1); 
            df3u2 = diff(vector_f{3}, u2); 
            df3u3 = diff(vector_f{3}, u3); 
            W = [double(subs((subs((subs(df1u1, u1, x0(1))), u2, x0(2))), u3, x0(3))) double(subs((subs((subs(df1u2, u1, x0(1))), u2, x0(2))), u3, x0(3))) double(subs((subs((subs(df1u3, u1, x0(1))), u2, x0(2))), u3, x0(3))); 
            double(subs((subs((subs(df2u1, u1, x0(1))), u2, x0(2))), u3, x0(3))) double(subs((subs((subs(df2u2, u1, x0(1))), u2, x0(2))), u3, x0(3))) double(subs((subs((subs(df2u3, u1, x0(1))), u2, x0(2))), u3, x0(3))); 
            double(subs((subs((subs(df3u1, u1, x0(1))), u2, x0(2))), u3, x0(3))) double(subs((subs((subs(df3u2, u1, x0(1))), u2, x0(2))), u3, x0(3))) double(subs((subs((subs(df3u3, u1, x0(1))), u2, x0(2))), u3, x0(3)))]; 
            W = inv(W); 
            v = W*f0; 
            x_next = x0 - v; 
        end 
        return; 
    else 
        x_next = [1000, 1000]; 
        while abs(x_next(1)-x0(1)) > 0.01 && abs(x_next(2)-x0(2)) > 0.01 
            c = c + 1; 
            if c > 5 
                break; 
            end 
            x_next = x0; 
            x0 = x_next; 
            f0 = [double(subs((subs((subs(vector_f{1}, u1, x0(1))), u2, x0(2))), t, t_now)); 
            double(subs((subs((subs(vector_f{2}, u1, x0(1))), u2, x0(2))), t, t_now))]; 
            df1u1 = diff(vector_f{1}, u1); 
            df1u2 = diff(vector_f{1}, u2); 
            df2u1 = diff(vector_f{2}, u1); 
            df2u2 = diff(vector_f{2}, u2); 
            W = [double(subs((subs((subs(df1u1, u1, x0(1))), u2, x0(2))), t, t_now)) double(subs((subs((subs(df1u2, u1, x0(1))), u2, x0(2))), t, t_now));  
            double(subs((subs((subs(df2u1, u1, x0(1))), u2, x0(2))), t, t_now)) double(subs((subs((subs(df2u2, u1, x0(1))), u2, x0(2))), t, t_now))]; 
            W = inv(W); 
            v = W*f0; 
            x_next = x0 - v; 
        end 
    end 
end

function explicit_method (count, u, tau_max, T, t, eps, f)
    y = u;
    if count == 2 
        f_1 = inline(f(1), 'u1', 'u2', 't'); 
        f_2 = inline(f(2), 'u1', 'u2', 't'); 
    else 
        f_1 = inline(f(1), 'u1', 'u2', 'u3', 't'); 
        f_2 = inline(f(2), 'u1', 'u2', 'u3', 't'); 
        f_3 = inline(f(3), 'u1', 'u2', 'u3', 't');
    end 
    while t < T 
        if count == 2 
            vector_f(1) = f_1(y(1), y(2), t); 
            vector_f(2) = f_2(y(1), y(2), t); 
        else 
            vector_f(1) = f_1(y(1), y(2), y(3), t); 
            vector_f(2) = f_2(y(1), y(2), y(3), t); 
            vector_f(3) = f_3(y(1), y(2), y(3), t); 
        end 
        for i = 1:count 
            tau(i) = eps/(abs(vector_f(i))+eps/tau_max); 
        end 
        step_tau = min(tau); 
        y = y + step_tau*vector_f; 
        hold on 
        plot(t, y(1), '.r', 'MarkerSize', 2); 
        plot(t, y(2), '.b', 'MarkerSize', 2); 
        if count == 3 
            plot(t, y(3), '.y', 'MarkerSize', 2); 
        end 
        t = t + step_tau; 
    end 
end

    
function implicit_method (count, u, tau_max, tau_min, T, eps, f) 
    t_now = tau_max; 
    t_next = t_now; 
    y_now = u; 
    y_back = u; 
    y_next = u; 
    tau_back = tau_min; 
    tau_now = tau_min; 
    f_1 = str2sym(f(1)); 
    f_2 = str2sym(f(2)); 
    if count == 3 
        f_3 = str2sym(f(3)); 
    end 
    while t_now < T 
        t_next = t_now + tau_now; 
        syms u1 u2 u3 t; 
        vector_f{1} = u1 - y_now (1) - tau_now*(f_1); 
        vector_f{2} = u2 - y_now (2) - tau_now*(f_2); 
        if count == 3 
            vector_f{3} = u3 - y_now (3) - tau_now*(f_3); 
        end 
        x0 = y_now'; 
        y_next = Newton(count, x0, vector_f, t_now)'; 
        eps_k = -(tau_now/(tau_now + tau_back))*(y_next - y_now - (tau_now/tau_back)*(y_now-y_back)); 
        for i = 1:count 
            if abs(eps_k(i)) > eps 
                tau_next_vector(i) = tau_now/2; 
            else 
                if eps/4 < abs(eps_k(i)) && abs(eps_k(i)) <= eps 
                    tau_next_vector(i) = tau_now; 
                else 
                    tau_next_vector(i) = 2*tau_now; 
                end 
            end 
        end 
        tau_next = min(tau_next_vector); 
        if tau_next > tau_max 
            tau_next = tau_max; 
        end 
        hold on 
        plot(t_now, y_next(1), '.r', 'MarkerSize', 2); 
        plot(t_now, y_next(2), '.b', 'MarkerSize', 2); 
        if count == 3 
            plot(t_now, y_next(3), '.black', 'MarkerSize', 2); 
        end 
    y_back = y_now; 
    y_now = y_next; 
    tau_back = tau_now; 
    tau_now = tau_next; 
    t_now = t_next; 
    end 
end 

function shihman_method (count, u, tau_max, tau, T, eps, f) 
    alpha_0 = 1; 
    alpha_1 = 0; 
    beta_0 = tau; 
    t_now = tau; 
    y_now = u; 
    y_back = u; 
    y_next = u; 
    tau_back = tau; 
    tau_now = tau; 
    f_1 = str2sym(f(1)); 
    f_2 = str2sym(f(2)); 
    if count == 3 
        f_3 = str2sym(f(3)); 
    end 
    c = 0; 
    while t_now < T 
        c = c + 1; 
        if c > 2 
            alpha_0 = ((tau_now+tau_back)^2)/(tau_back*(2*tau_now+tau_back)); 
            alpha_1 = -((tau_now)^2)/(tau_back*(2*tau_now+tau_back));
            beta_0 = tau_now*(tau_now+tau_back)/(2*tau_now+tau_back); 
        end 
        t_next = t_now + tau_now; 
        syms u1 u2 u3 t;
        vector_f{1} = u1 - alpha_1 * y_back (1) - alpha_0 * y_now (1) - beta_0 * (f_1);
        vector_f{2} = u2 - alpha_1 * y_back (2) - alpha_0 * y_now (2) - beta_0 *(f_2); 
        if count == 3
            vector_f{3} = u3 - alpha_1 * y_back (3) - alpha_0 * y_now (3) - beta_0 * (f_3); 
        end 
        x0 = y_now'; 
        y_next = Newton(count, x0, vector_f, t_now)';
        eps_k = -(tau_now/(tau_now + tau_back))*(y_next - y_now - 192 (tau_now/tau_back)*(y_now-y_back)); 
        for i = 1:count 
            if abs(eps_k(i)) > eps 
                tau_next_vector(i) = tau_now/2; 
            else 
                if eps/4 < abs(eps_k(i)) && abs(eps_k(i)) <= eps 
                    tau_next_vector(i) = tau_now; 
                else 
                    tau_next_vector(i) = 2*tau_now; 
                end 
            end 
        end 
        tau_next = min(tau_next_vector); 
        if tau_next > tau_max 
            tau_next = tau_max; 
        end 
        hold on 
        plot(t_now, y_next(1), '.r', 'MarkerSize', 2); 
        plot(t_now, y_next(2), '.b', 'MarkerSize', 2); 
        if count == 3 
            plot(t_now, y_next(3), '.black', 'MarkerSize', 2); 
        end 
        y_back = y_now; 
        y_now = y_next; 
        tau_back = tau_now; 
        tau_now = tau_next; 
        t_now = t_next; 
    end 
end 