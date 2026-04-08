clear 
clc

%% parameters
c = 1.5;            % root chord
lambda = 2/3;       % taper
L = sqrt(2);        % semi span
S = 2*sqrt(2);      % planform area
Vinf = 50;          % airspeed m/s
rho = 1.225;        % density kg/m^3
adeg = 6;           % root angle of attack degrees
alpha = adeg*pi/180; 
NP = 80;            % spanwise panels (reduced slightly for speed)

% sweep angles 
sweeps = [45, 0, -45];
labels = {'45^\circ, Aft (Numerical)', '0^\circ, Unswept (Numerical)', '-45^\circ, Forward (Numerical)'};

% initialize figures
fig1 = figure('Name', 'Cl Distribution'); hold on; grid on;
fig2 = figure('Name', 'Spanwise Loading'); hold on; grid on;

%% main loop over sweep angles
for s = 1:length(sweeps)
    LAM = sweeps(s);

    % mesh
    DY = 2/NP;
    Y_L = 0:DY:1;
    y_L = -1:DY:1;
    C_Y = c*(1-lambda*Y_L);
    c_Y = [C_Y(end:-1:2) C_Y];
    
    x_hc = sign(y_L).*y_L*L*tand(LAM);    % half chord
    x_qc = -0.25*c_Y + x_hc;              % quarter chord
    x_3qc = 0.25*c_Y + x_hc;              % 3/4 chord
    z = zeros(size(x_hc));
    
    % collocation points
    A = [x_qc(1:NP); y_L(1:NP)*L; z(1:NP)];        % Left vortex corner
    B = [x_qc(2:NP+1); y_L(2:NP+1)*L; z(2:NP+1)];  % Right vortex corner
    Cy = 0.5*(A(2,:) + B(2,:));
    Cx = interp1(y_L, x_3qc, Cy/L);
    C = [Cx; Cy; zeros(1, NP)];
    
    AIC_mat = zeros(NP, NP);
    for j = 1:NP
        for k = 1:NP
            VAB = V_AB(A(:,k), B(:,k), C(:,j));
            VAI = VA_INF(A(:,k), C(:,j));
            VBI = VB_INF(B(:,k), C(:,j));
            panel_n = cross(A(:,k)-C(:,k), A(:,k)-B(:,k));
            panel_n = panel_n / norm(panel_n);            
            V_total = VAB + VAI + VBI;
            AIC_mat(j,k) = dot(V_total, panel_n);
        end
    end
    
    % normal vector for panels 
    n_global = [0; 0; 1];               % approximate as vertical 
    V_rhs = -Vinf * sin(alpha); 
    gamma = AIC_mat \ (V_rhs * ones(NP, 1));
    
    % lift
    qinf = 0.5 * rho * Vinf^2;
    Li = rho * Vinf * gamma' .* (DY * L);
    Ltot = sum(Li);
    Cltot = Ltot / (qinf * S);

    % one side of wing 
    Li_half = Li((NP/2 + 1):NP);
    
    % local chord and local Cl
    mid_y = Cy((NP/2 + 1):NP);
    local_c = 0.5 * (c_Y((NP/2 + 1):NP) + c_Y((NP/2 + 2):NP+1));
    Cl_local = Li_half ./ (qinf * local_c * (DY * L));
    c_avg = S/(2*L);

    % plotting 1: CL 
    figure(fig1);
    plot(mid_y/L, Cl_local/Cltot, '-', 'LineWidth', 2.5, 'DisplayName', labels{s});
    
    % plotting 2: loading ---
    figure(fig2);
    plot(mid_y/L, (Cl_local.*local_c)/(c_avg*Cltot), '-', 'LineWidth', 2.5, 'DisplayName', labels{s});
end

%% reference data
% aft sweep Cl
cxa = [0.017; 0.064; 0.108; 0.152; 0.200; 0.244; 0.291; 0.335; 0.379; 0.427; ...
     0.471; 0.517; 0.563; 0.609; 0.653; 0.698; 0.745; 0.788; 0.836; 0.882; ...
     0.927; 0.975];
cya = [0.799; 0.830; 0.857; 0.888; 0.915; 0.947; 0.974; 0.999; 1.032; 1.055; ...
     1.086; 1.110; 1.135; 1.154; 1.172; 1.187; 1.189; 1.185; 1.157; 1.095; ...
     0.972; 0.726];
% unswept Cl
cx0 = [0.020; 0.063; 0.108; 0.155; 0.197; 0.246; 0.291; 0.335; 0.384; 0.426; ...
     0.474; 0.519; 0.566; 0.609; 0.654; 0.701; 0.745; 0.793; 0.838; 0.882; ...
     0.931; 0.976];
cy0 = [0.873; 0.902; 0.924; 0.941; 0.968; 0.987; 1.009; 1.026; 1.046; 1.067; ...
     1.075; 1.096; 1.102; 1.115; 1.115; 1.115; 1.101; 1.081; 1.036; 0.962; ...
     0.847; 0.624];
% fore sweep Cl
cxf = [0.019; 0.063; 0.107; 0.155; 0.199; 0.245; 0.290; 0.335; 0.384; 0.426; ...
     0.474; 0.519; 0.566; 0.610; 0.653; 0.702; 0.745; 0.794; 0.841; 0.883; ...
     0.931; 0.977];
cyf = [0.961; 0.988; 1.004; 1.025; 1.040; 1.050; 1.065; 1.071; 1.069; 1.075; ...
     1.071; 1.063; 1.057; 1.043; 1.029; 1.007; 0.978; 0.939; 0.890; 0.820; ...
     0.700; 0.515];

% aft sweep L
lxa = [0.026, 0.071, 0.116, 0.161, 0.206, 0.252, 0.296, 0.341, 0.389, 0.433, ...
    0.48, 0.524, 0.572, 0.616, 0.66, 0.706, 0.751, 0.799, 0.843, 0.887, 0.934, 0.978];
lya = [1.181, 1.189, 1.192, 1.187, 1.187, 1.179, 1.175, 1.162, 1.149, 1.129,...
    1.106, 1.081, 1.058, 1.028, 0.99, 0.955, 0.899, 0.84, 0.767, 0.68, 0.561, 0.382];
% unswept L
lx0 = [0.026, 0.071, 0.113, 0.161, 0.206, 0.254, 0.299, 0.342, 0.389, 0.434, 0.480, 0.524,...
    0.572, 0.614, 0.660, 0.707, 0.751, 0.799, 0.842, 0.887, 0.933, 0.979];
ly0 = [1.291, 1.293, 1.276, 1.271, 1.256, 1.238, 1.218, 1.193, 1.165, 1.134, 1.106, 1.066,...
    1.030, 0.987, 0.944, 0.891, 0.832, 0.764, 0.693, 0.596, 0.487, 0.329];
% fore sweep L
lxf = [0.026, 0.070, 0.114, 0.162, 0.206, 0.252, 0.297, 0.342, 0.388, 0.434, 0.481, 0.525, 0.572, ...
    0.615, 0.660, 0.707, 0.750, 0.799, 0.843, 0.888, 0.933, 0.979];
lyf = [1.420, 1.408, 1.393, 1.380, 1.347, 1.320, 1.282, 1.236, 1.200, 1.150, 1.099, 1.048, 0.990, ...
    0.929, 0.865, 0.802, 0.738, 0.672, 0.591, 0.510, 0.408, 0.273];

%% Plot 1
figure(fig1)
hold on 
plot(cxa, cya, '^k', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'Aft (Reference)')
plot(cx0, cy0, 'sk', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Unswept (Reference)')
plot(cxf, cyf, 'xk', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'Forward (Reference)')
title('Effect of Sweep on Span-wise Lift Distribution', 'FontSize', 30);
xlabel('Non-dimensional Spanwise Position (y/L)', 'FontSize', 25);
ylabel('Normalised Lift Coefficient, C_l / C_{L,Total}', 'FontSize', 25);
legend('Location', 'southwest', 'FontSize', 20);
ylim([0.4 1.3]);
set(gca, 'fontsize', 20)
hold off;

%% Plot 2 
figure(fig2);
hold on 
plot(cxa, lya, '^k', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'Aft (Reference)')
plot(cx0, ly0, 'sk', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Unswept (Reference)')
plot(cxf, lyf, 'xk', 'MarkerSize', 15, 'LineWidth', 1.5, 'DisplayName', 'Forward (Reference)')
title('Spanwise Loading Distribution', 'FontSize', 30);
xlabel('Non-dimensional Spanwise Position (y/L)', 'FontSize', 25);
ylabel('Normalised Loading, c\cdotc_l/c_a', 'FontSize', 25);
legend('Location', 'southwest', 'FontSize', 20);
grid on; set(gca, 'fontsize', 20);