% Simple VLM for general planforms, comparing multiple sweep angles
clear; clc;

%% Parameters & Constants
c = 1.5;            % root chord
lambda = 2/3;       % taper
L = sqrt(2);        % semi span
S = 2*sqrt(2);      % planform area
Vinf = 50;          % airspeed m/s
rho = 1.225;        % density kg/m^3
adeg = 6;           % root angle of attack degrees
alpha = adeg*pi/180; 
NP = 60;            % spanwise panels (reduced slightly for speed)

% Sweep angles to iterate through
sweeps = [45, 0, -45];
labels = {'45^\circ, Aft (Numerical)', '0^\circ, Unswept (Numerical)', '-45^\circ, Forward (Numerical)'};

% Initialize Figure for Cl Distribution
figure('Name', 'Sweep Comparison');
hold on; grid on;

%% Main Loop over Sweep Angles
for s = 1:length(sweeps)
    LAM = sweeps(s);
    
    % --- Mesh Generation ---
    DY = 2/NP;
    Y_L = 0:DY:1;
    y_L = -1:DY:1;
    C_Y = c*(1-lambda*Y_L);
    c_Y = [C_Y(end:-1:2) C_Y];
    
    x_hc = sign(y_L).*y_L*L*tand(LAM);    % half chord
    x_qc = -0.25*c_Y + x_hc;              % quarter chord
    x_3qc = 0.25*c_Y + x_hc;              % 3/4 chord
    z = zeros(size(x_hc));
    
    % --- Vortex Points & Collocation Points ---
    A = [x_qc(1:NP); y_L(1:NP)*L; z(1:NP)];        % Left vortex corner
    B = [x_qc(2:NP+1); y_L(2:NP+1)*L; z(2:NP+1)];  % Right vortex corner
    Cy = 0.5*(A(2,:) + B(2,:));
    Cx = interp1(y_L, x_3qc, Cy/L);
    C = [Cx; Cy; zeros(1, NP)];
    
    % --- Aerodynamic Influence Coefficients (AIC) ---
    AIC_mat = zeros(NP, NP);
    for j = 1:NP
        for k = 1:NP
            % These functions (V_AB, VA_INF, VB_INF) must be in your path/folder
            VAB = V_AB(A(:,k), B(:,k), C(:,j));
            VAI = VA_INF(A(:,k), C(:,j));
            VBI = VB_INF(B(:,k), C(:,j));
            
            % Unit normal n calculation (simplified for flat plate)
            % Assuming n = [sin(alpha), 0, cos(alpha)] is too simple if wing is complex,
            % but for this case, we use the panel normal:
            panel_n = cross(A(:,k)-C(:,k), A(:,k)-B(:,k));
            panel_n = panel_n / norm(panel_n);
            
            V_total = VAB + VAI + VBI;
            AIC_mat(j,k) = dot(V_total, panel_n);
        end
    end
    
    % --- Solve System ---
    % Boundary condition: V_inf dot n
    % Normal vector for panels (approximate as vertical for small alpha)
    n_global = [0; 0; 1]; 
    V_rhs = -Vinf * sin(alpha); 
    gamma = AIC_mat \ (V_rhs * ones(NP, 1));
    
    % --- Lift Calculation ---
    qinf = 0.5 * rho * Vinf^2;
    Li = rho * Vinf * gamma' .* (DY * L);
    Ltot = sum(Li);
    Cltot = Ltot / (qinf * S);

    % data for one side of wing (right)
    Li_half = Li((NP/2 + 1):NP);
    
    % Local Chord and local Cl
    mid_y = Cy((NP/2 + 1):NP);
    local_c = 0.5 * (c_Y((NP/2 + 1):NP) + c_Y((NP/2 + 2):NP+1));
    Cl_local = Li_half ./ (qinf * local_c * (DY * L));
    
    % --- Plotting Results ---
    plot(mid_y/L, Cl_local/Cltot, '-', 'LineWidth', 2.5, 'DisplayName', labels{s});
end

%% reference data
% aft sweep
xa = [0.017; 0.064; 0.108; 0.152; 0.200; 0.244; 0.291; 0.335; 0.379; 0.427; ...
     0.471; 0.517; 0.563; 0.609; 0.653; 0.698; 0.745; 0.788; 0.836; 0.882; ...
     0.927; 0.975];
ya = [0.799; 0.830; 0.857; 0.888; 0.915; 0.947; 0.974; 0.999; 1.032; 1.055; ...
     1.086; 1.110; 1.135; 1.154; 1.172; 1.187; 1.189; 1.185; 1.157; 1.095; ...
     0.972; 0.726];

% unswept
x0 = [0.020; 0.063; 0.108; 0.155; 0.197; 0.246; 0.291; 0.335; 0.384; 0.426; ...
     0.474; 0.519; 0.566; 0.609; 0.654; 0.701; 0.745; 0.793; 0.838; 0.882; ...
     0.931; 0.976];
% unswept
y0 = [0.873; 0.902; 0.924; 0.941; 0.968; 0.987; 1.009; 1.026; 1.046; 1.067; ...
     1.075; 1.096; 1.102; 1.115; 1.115; 1.115; 1.101; 1.081; 1.036; 0.962; ...
     0.847; 0.624];

% fore sweep
xf = [0.019; 0.063; 0.107; 0.155; 0.199; 0.245; 0.290; 0.335; 0.384; 0.426; ...
     0.474; 0.519; 0.566; 0.610; 0.653; 0.702; 0.745; 0.794; 0.841; 0.883; ...
     0.931; 0.977];
yf = [0.961; 0.988; 1.004; 1.025; 1.040; 1.050; 1.065; 1.071; 1.069; 1.075; ...
     1.071; 1.063; 1.057; 1.043; 1.029; 1.007; 0.978; 0.939; 0.890; 0.820; ...
     0.700; 0.515];

%% Finalize Plot
hold on 
plot(xa, ya, '^k', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'Aft (Reference)')
plot(x0, y0, 'sk', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Unswept (Reference)')
plot(xf, yf, 'xk', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'Forward (Reference)')
title('Effect of Sweep on Span-wise Lift Distribution', 'FontSize', 30);
xlabel('Non-dimensional Spanwise Position (y/L)', 'FontSize', 25);
ylabel('C_l / C_{L,Total}', 'FontSize', 25);
legend('Location', 'southwest', 'FontSize', 20);
ylim([0.4 1.3]);
set(gca, 'fontsize', 20)
hold off;
