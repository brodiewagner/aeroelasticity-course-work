function [Li, Cltot] = solve_VLM(alpha_vector, Vinf, rho, S, NP, A, B, C, n, DY, L)
    
     for j = 1:NP
        for k = 1:NP
            VAB = V_AB(A(:,k),B(:,k),C(:,j)) ;
            VAI = VA_INF(A(:,k),C(:,j)) ;
            VBI = VB_INF(B(:,k),C(:,j)) ;
            AIC(j,k,:) = VAB+VAI+VBI ;
        end
     end
     
    AICx = squeeze(AIC(:,:,1)) ;    % get x component of downwash
    AICy = squeeze(AIC(:,:,2)) ;    % get y component of downwash
    AICz = squeeze(AIC(:,:,3)) ;    % get z component of downwash
    
    % form aerodynamic influence coefficient matrix projection of V on n
    AIC_mat = zeros(NP, NP) ;
    for i = 1:NP
        for j = 1:NP
            AIC_mat(i,j) = n(1,j)*AICx(i,j) + n(2,j)*AICy(i,j) + n(3,j)*AICz(i,j) ;     % this is just a dot product
        end
    end
    
    % 2. Calculate local incident flow for each panel
    % Use the specific alpha for each panel from the elastic solver
    V = - Vinf*(cos(alpha_vector).*n(1,:) + sin(alpha_vector).*n(3,:));
    
    % 3. Solve for Gamma
    gamma = AIC_mat \ V';
    
    % 4. Compute Local Lift (Li)
    Li = rho*Vinf*gamma'.*(DY*L) ;   
    Li = Li(:, 1) ;                       % force 1×NP row vector
    
    % 5. Total Lift Coefficient
    Ltot = sum(Li);
    qinf = 0.5 * rho * Vinf^2;
    Cltot = Ltot / (qinf * S);
end