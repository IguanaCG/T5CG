clear all
close all
clc

% PARAMETROS 
K = [5, 0, 0; 0, 5, 0; 0, 0, 5];
dt = 0.05;
S = 30;
t = 0 : dt : S;
q0 = [0.1; 0.2; 0.1; 0.0];

% TRAYECTORIAS 
funcion_trayectoria = @(ti) trayectoria_circulo(ti, S);
%funcion_trayectoria = @(ti) trayectoria_rosa(ti, S);
%funcion_trayectoria = @(ti) trayectoria_triangulo(ti, S);

% SIMULACION 
[hist_t, hist_x, hist_xd, hist_q, hist_qdot] = ejecutarControlador(t, dt, q0, K, funcion_trayectoria);
graficarResultados(hist_t, hist_x, hist_xd, hist_q, hist_qdot);

% FUNCIONES 
function [hist_t, hist_x, hist_xd, hist_q, hist_qdot] = ejecutarControlador(t, dt, q, K, f_tray)
    N = length(t);
    hist_t = zeros(1, N); hist_x = zeros(3, N); hist_xd = zeros(3, N);
    hist_q = zeros(4, N); hist_qdot = zeros(4, N);
    
    figure(1); hold on; grid on; axis([-0.6 0.8 -0.6 0.8 0 1.2]);
    view([0 90]);
    
    contador = 0;
    for i = 1:N
        tiempo = t(i);
        [xd, vxd] = f_tray(tiempo);
        
        [A1, A2, A3, A4] = obtenerMatricesDH(q);
        T04 = A1 * A2 * A3 * A4;
        x_act = T04(1:3, 4);
        
        u = vxd + K * (xd - x_act);
        
        % Llamada al Jacobiano Analítico
        J = calcularJacobiano(q, x_act);
        q_dot = pinv(J) * u;
        q = q + q_dot * dt;
        
        hist_t(i) = tiempo; hist_x(:,i) = x_act; hist_xd(:,i) = xd;
        hist_q(:,i) = q; hist_qdot(:,i) = q_dot;
        
        contador = contador + 1;
        if contador == 5
            cla;
            plot3(hist_xd(1,1:i), hist_xd(2,1:i), hist_xd(3,1:i), 'r', 'LineWidth', 2);
            Dibujar_Manipulador({A1, A2, A3, A4}, {'RRPR'}, 0.1);
            view([0 90]); axis([-0.6 0.8 -0.6 0.8 0 1.2]); 
            drawnow; contador = 0;
        end
    end
end

% JACOBIANO 
function J = calcularJacobiano(q, x_act)
    q1 = q(1);
    q2 = q(2);
   
    J11 = -0.15 * sin(q1) - 0.4 * sin(q1 + q2);
    J12 = -0.4 * sin(q1 + q2);
    J21 = 0.15 * cos(q1) + 0.4 * cos(q1 + q2);
    J22 = 0.4 * cos(q1 + q2);
    J33 = -1;

    J = [ J11,  J12,   0,   0;
          J21,  J22,   0,   0;
            0,    0, J33,   0 ];
end

function [A1, A2, A3, A4] = obtenerMatricesDH(q)
    A1 = [cos(q(1)) -sin(q(1)) 0 0.15*cos(q(1)); sin(q(1)) cos(q(1)) 0 0.15*sin(q(1)); 0 0 1 0.8; 0 0 0 1];
    A2 = [cos(q(2)) sin(q(2)) 0 0.4*cos(q(2)); sin(q(2)) -cos(q(2)) 0 0.4*sin(q(2)); 0 0 -1 0; 0 0 0 1];
    A3 = [1 0 0 0; 0 1 0 0; 0 0 1 q(3)+0.3; 0 0 0 1];
    A4 = [cos(q(4)) -sin(q(4)) 0 0; sin(q(4)) cos(q(4)) 0 0; 0 0 1 0.2; 0 0 0 1];
end

% TRAYECTORIAS 
function [xd, vxd] = trayectoria_circulo(tiempo, S)
    r = 0.1; cX = 0.2; cY = 0.35; cZ = 0.6;
    w = (2*pi)/S;
    xd = [cX + r*cos(w*tiempo); cY + r*sin(w*tiempo); cZ];
    vxd = [-r*w*sin(w*tiempo); r*w*cos(w*tiempo); 0];
end

function [xd, vxd] = trayectoria_rosa(tiempo, S)
    a = 0.07; k = 8; cX = 0.3; cY = 0.35; cZ = 0.6;
    w = (2*pi)/S; th = w*tiempo;
    r = a*cos(k*th);
    dr = -a*k*w*sin(k*th);
    xd = [cX + r*cos(th); cY + r*sin(th); cZ];
    vxd = [dr*cos(th) - r*w*sin(th); dr*sin(th) + r*w*cos(th); 0];
end

function [xd, vxd] = trayectoria_triangulo(tiempo, S)
    P1 = [0.30; 0.25; 0.6]; 
    P2 = [0.35; 0.15; 0.6]; 
    P3 = [0.25; 0.15; 0.6]; 
    T_l = S/3;
    if tiempo < T_l
        Pi = P1; Pf = P2; t_loc = tiempo;
    elseif tiempo < 2*T_l
        Pi = P2; Pf = P3; t_loc = tiempo - T_l;
    else
        Pi = P3; Pf = P1; t_loc = tiempo - 2*T_l;
    end
    xd = Pi + (t_loc/T_l)*(Pf - Pi);
    vxd = (1/T_l)*(Pf - Pi);
end

% --- GRAFICAS ---
function graficarResultados(t, x_act, xd, q, qdot)
    figure(2); 
    clf; 
    hold on; 
    grid on;
    
    % Eje X 
    plot(t, xd(1,:), '*y', 'LineWidth', 0.3); % X Deseada
    plot(t, x_act(1,:), '-b', 'LineWidth', 1.5); % X Actual
    
    % Eje Y 
    plot(t, xd(2,:), '*y', 'LineWidth', 0.3); % Y Deseada
    plot(t, x_act(2,:), '-g', 'LineWidth', 1.5); % Y Actual
    
    % Eje Z 
    plot(t, xd(3,:), '*y', 'LineWidth', 0.3); % Z Deseada
    plot(t, x_act(3,:), '-r', 'LineWidth', 1.5); % Z Actual
    
    title('Posición Deseada vs Actual (Ejes X, Y, Z)');
    xlabel('Tiempo (s)');
    ylabel('Posición (m)');
    
    legend('X Deseada', 'X Actual', 'Y Deseada', 'Y Actual', 'Z Deseada', 'Z Actual', 'Location', 'best');
    
    % 2. Gráfica de Velocidades Articulares
    figure(3); 
    plot(t, qdot, 'LineWidth', 1.5); 
    title('Velocidades Articulares'); 
    xlabel('Tiempo (s)'); ylabel('Velocidad');
    legend('q1','q2','d3','q4'); 
    grid on;
    
    figure(4); 
    plot(t, q, 'LineWidth', 1.5); 
    title('Posiciones Articulares'); 
    xlabel('Tiempo (s)'); ylabel('Posición');
    legend('q1','q2','d3','q4'); 
    grid on;
end