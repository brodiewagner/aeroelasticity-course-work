function [Li, Cltot] = solve_VLM_2(alpha_row, Vinf, rho, S, NP, A, B, C, n, DY, L)
% solve_VLM  Vortex Lattice Method solver
%
% INPUTS:
%   alpha_row  - 1xNP row vector of local angles of attack (rad)
%   Vinf       - freestream speed (m/s)
%   rho        - air density (kg/m^3)
%   S          - reference wing area (m^2)
%   NP         - number of panels
%   A, B       - 3xNP vortex corner points (quarter chord)
%   C          - 3xNP control points (3/4 chord)
%   n          - 3xNP unit normals at control points
%   DY         - dimensionless panel span
%   L          - semi-span (m)
%
% OUTPUTS:
%   Li     - 1xNP lift per panel (N)
%   Cltot  - total lift coefficient

    % Force alpha to be a 1xNP row vector
    alpha_row = reshape(alpha_row, 1, NP) ;

    % --- Build AIC matrix ---
    % AIC(j,k) = normal component of velocity induced at control point j
    %            by a unit-strength horseshoe vortex on panel k
    AIC_mat = zeros(NP, NP) ;

    for j = 1:NP          % control point (where we apply BC)
        for k = 1:NP      % vortex panel (source of induced velocity)
            vel = V_AB(A(:,k), B(:,k), C(:,j)) ...
                + VA_INF(A(:,k), C(:,j))       ...
                + VB_INF(B(:,k), C(:,j)) ;     % 3x1 velocity vector

            % Dot with normal at control point j
            AIC_mat(j,k) = dot(n(:,j), vel) ;
        end
    end

    % --- RHS: normal component of freestream at each control point ---
    % V_inf direction: [cos(alpha), 0, sin(alpha)] for each panel
    % Normal BC: (V_inf + V_induced) . n = 0  =>  V_induced.n = -V_inf.n
    RHS = zeros(NP, 1) ;
    for j = 1:NP
        V_inf_j = Vinf * [cos(alpha_row(j)) ; 0 ; sin(alpha_row(j))] ;
        RHS(j)  = -dot(n(:,j), V_inf_j) ;
    end

    % --- Solve for circulation strengths ---
    gamma = AIC_mat \ RHS ;     % NP x 1

    % --- Lift per panel (Kutta-Joukowski, one wing) ---
    Li = rho * Vinf * gamma' * (DY*L) ;    % 1xNP row vector

    % --- Total lift coefficient ---
    Ltot  = sum(Li) ;
    qinf  = 0.5*rho*Vinf^2 ;
    Cltot = Ltot / (qinf*S) ;

end