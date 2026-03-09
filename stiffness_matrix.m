function[E] = stiffness_matrix(N,M,EI,GJ,L,y_L)

% for ii=1:N
%     psi(ii,:) =(y_L).^(ii+1); % ith bending function wing
%     psid(ii,:)=(ii+1)*(y_L).^(ii);
% end
% for ii=1:M
%     phi(ii,:)= (y_L).^(ii); % ith torsion function wing
% end

  for i = 1:N
        psii_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second
        for j=1:N
           psij_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;        % jth wing bending function second derivative
           % Del, Dels and B matrices
            B(i,j) = (L-0)*trapz(y_L,psii_double.*psij_double);
        end
    end
    % setup i and j torsion for wing and store,2 torsion modes
    for i = 1:M
            phii_single = (i/L).*((y_L).^(i-1)) ;            % ith torsion function first derivative wing
            for j = 1:M
                 phij_single = (j/L).*((y_L).^(j-1)) ; % j-th torsion function first derivative wing
                 % determine the D, D_store and T matrices
                 T(i,j) = (L-0)*trapz(y_L,phii_single.*phij_single) ;
            end
    end
 
    % matrix with no coupling between torsion and bending
    E=[EI*B,         zeros(N,M) ;
       (zeros(N,M))',  GJ*T   ] ;

end