clear ;
clc ;

%% panels
NP_s_list = [10, 30, 70, 100, 200] ;   % half-span structural panels
n_cases   = numel(NP_s_list) ;

%% parameters
c = 1.5 ;         % root chord (m)
lambda = 0.4 ;    % taper ratio
L = 6 ;           % semi-span (m)
e = 0.25 ;        % elastic axis offset (chords aft of quarter-chord)
LAM = 0 ;        % half-chord sweep angle (deg)

mg = 10e3*9.81 ;
Vinf = 150 ;
rho = 0.5238*1.225 ;
qinf = 0.5*rho*Vinf^2 ;

EI = 2e6 ;
GJ = 5e5 ;
N  = 3   ;   % bending modes
M  = 2   ;   % torsion modes

%% storage
y_common   = linspace(0, 1, 201) ;
CL_store   = zeros(n_cases, numel(y_common)) ;
alpha_trim = zeros(1, n_cases) ;

%% loop
for ic = 1:n_cases
    NP_s = NP_s_list(ic) ;
    fprintf('\n=== NP_s = %d ===\n', NP_s) ;

    % structural mesh 
    DY_s = 1/NP_s ;
    y_L  = 0:DY_s:1 ;
    c_Y  = c*(1 - lambda*y_L) ;
    dy   = DY_s*L ;
    S    = trapz(c_Y)*dy*2 ;
    AR   = (2*L)^2/S ;

    % VLM mesh (full-span)
    NP_v  = 2*NP_s ;
    DY_v  = 2/NP_v ;
    Y_L_v = 0:DY_v:1 ;
    y_L_v = -1:DY_v:1 ;

    C_Y_v = c*(1 - lambda*Y_L_v) ;
    c_Y_v = [C_Y_v(end:-1:2), C_Y_v] ;

    x_hc_v = abs(y_L_v) * L * tand(LAM) ;
    x_qc_v  = -0.25*c_Y_v + x_hc_v ;
    x_3qc_v =  0.25*c_Y_v + x_hc_v ;
    z_v     = zeros(size(y_L_v)) ;

    A_v = [x_qc_v(1:NP_v);     y_L_v(1:NP_v)*L;     z_v(1:NP_v)    ] ;
    B_v = [x_qc_v(2:NP_v+1);   y_L_v(2:NP_v+1)*L;   z_v(2:NP_v+1)  ] ;

    Cy_v = 0.5*(A_v(2,:) + B_v(2,:)) ;
    Cx_v = interp1(y_L_v, x_3qc_v, Cy_v/L) ;
    C_v  = [Cx_v; Cy_v; 0*Cy_v] ;

    n_v = zeros(3, NP_v) ;
    for k = 1:NP_v
        nk = cross(A_v(:,k) - C_v(:,k), A_v(:,k) - B_v(:,k)) ;
        n_v(:,k) = nk / norm(nk) ;
    end

    E_stiff = stiffness_matrix(N, M, EI, GJ, L, y_L) ;

    psi_i  = zeros(N, NP_s+1) ;
    psi_id = zeros(N, NP_s+1) ;
    phi_i  = zeros(M, NP_s+1) ;
    for ii = 1:N
        psi_i(ii,:)  = y_L.^(ii+1) ;
        psi_id(ii,:) = ((ii+1)/L)*y_L.^ii ;
    end
    for ii = 1:M
        phi_i(ii,:) = y_L.^ii ;
    end

    a3 = 2*pi*AR / (2 + sqrt(4 + AR^2)) ;
    alpha_cur = mg/(qinf*S*a3) ;

    Li_s  = qinf*c_Y*a3*alpha_cur ;
    Mi_s  = e*c*Li_s ;

    F = zeros(N+M, 1) ;
    for ii = 1:N;  F(ii)   = trapz(y_L, Li_s.*psi_i(ii,:))*L;  end
    for ii = 1:M;  F(ii+N) = trapz(y_L, Mi_s.*phi_i(ii,:))*L;  end

    eta   = E_stiff\F ;
    theta = phi_i' * eta(N+1:N+M) ;
    wd    = psi_id'* eta(1:N) ;

    y_panels = 0.5*(y_L(1:end-1) + y_L(2:end)) ;

    %% Trim loop
    for trimstep = 2:300

        alphae_half   = alpha_cur + theta*cosd(LAM) - wd*sind(LAM) ;
        alphae_panels = 0.5*(alphae_half(1:end-1) + alphae_half(2:end)) ;
        alphae_full   = [flip(alphae_panels); alphae_panels]' ;

        [Li_vlm] = solve_VLM(alphae_full, Vinf, rho, S, NP_v, A_v, B_v, C_v, n_v, DY_v, L) ;

        Li_star = Li_vlm(NP_v/2+1 : NP_v) / dy ;
        Li_s    = interp1(y_panels, Li_star, y_L, 'linear', 'extrap') ;

        Ltot  = sum(Li_vlm) ;
        Lover = Ltot - mg ;
        aover = Lover / (qinf*S*a3) ;
        alpha_cur = alpha_cur - aover*0.01 ;

        Mi_s = e*c*Li_s ;
        F = zeros(N+M, 1) ;
        for ii = 1:N;  F(ii)   = trapz(y_L, Li_s.*psi_i(ii,:))*L;  end
        for ii = 1:M;  F(ii+N) = trapz(y_L, Mi_s.*phi_i(ii,:))*L;  end

        eta   = E_stiff\F ;
        theta = phi_i'*eta(N+1:N+M) ;
        wd    = psi_id'*eta(1:N) ;

        if abs(Lover) < mg*0.0005
            fprintf('  Converged at step %d (alpha=%.5f rad)\n', trimstep, alpha_cur) ;
            break
        end
    end

    % store result 
    CL_y = Li_s ./ (qinf*c_Y) ;
    CL_store(ic,:) = interp1(y_L, CL_y, y_common, 'linear', 'extrap') ;
    alpha_trim(ic) = alpha_cur ;

end 

%% plot 1
colors = lines(n_cases) ;
figure('Name','CL_y – Panel Convergence Study', 'Position',[100 100 900 520])
hold on; grid on; box on
for ic = 1:n_cases
    plot(y_common, CL_store(ic,:), ...
         'Color', colors(ic,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('N_{ps} = %d', NP_s_list(ic))) ;
end
xlabel('Dimensionless span  y/L',       'FontSize', 14)
ylabel('Local lift coefficient  C_L',   'FontSize', 14)
title('Spanwise C_L Convergence with Panel Resolution', 'FontSize', 14)
legend('Location','southwest', 'FontSize', 11)

%% plot 2
stations     = [0.25, 0.50, 0.75, 0.95] ;
station_lbls = {'y/L = 0.25','y/L = 0.50','y/L = 0.75','y/L = 0.95'} ;

CL_at_stations = zeros(n_cases, numel(stations)) ;
for is = 1:numel(stations)
    [~,idx] = min(abs(y_common - stations(is))) ;
    CL_at_stations(:,is) = CL_store(:,idx) ;
end

figure('Name','CL Point Convergence', 'Position',[120 120 900 520])
hold on; grid on; grid minor; box on
station_colors = lines(numel(stations)) ;
for is = 1:numel(stations)
    plot(NP_s_list, CL_at_stations(:,is), '-o', ...
         'Color', station_colors(is,:), 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', station_lbls{is}) ;
end
set(gca,'XTick', NP_s_list, 'XTickLabel', string(NP_s_list))
xlim([10, 200])
xlabel('Number of panels  N_{p}', 'FontSize', 14)
ylabel('Local lift coefficient  C_L', 'FontSize', 14)
title('C_L Point Convergence at Spanwise Stations', 'FontSize', 14)
legend('Location','east', 'FontSize', 11)

