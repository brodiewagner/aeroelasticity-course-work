clear 
clc

% wing parameters
mw=;
c=;
b=;
L=;
EA=;
IA=;
x_b=;
EI=;
GJ=;
mbar=;
Iwing=;
rho=mbar/(pi*b*b*32.6);
e=(EA/c)-(1/4);
a=(e*c)/b-0.5; % elastic axis position

% parameters for wing store 4 
ms = 0.636*mw ;           % store mass, kg
xs = -0.687 ;             % position of store ahead of Elastic Axis, m
Is = 2.68*Iwing*L ;       % pitch inertia of store about centre , kg m^2
x = xs*b ;

% store position at 100 points along span
for int_pos = 1:1:101
    Ls = (int_pos-1).*(L/100) ;
    y_L = 0:0.01:1 ;
    N = 3 ;
    M = 2 ;
    % calculate Delta, Deltas, B, C, D and T matrices setup functions for bending, 3 bending and 2 torsion modes 
   for i=1:N
        psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
        psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
        for j = 1:N
           psi_j = (y_L).^(j+1) ;                                   % jth bedning function of wing
           psi_j_store = (Ls/L).^(j+1) ;                            % jth bending mode of store
           psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;        % jth wing bending function second derivative (psi'') 
           % Del, Dels and B matrices
           Del(i,j) = (L-0)*trapz(y_L,psi_i.*psi_j) ;              % Delta for wing
           Dels(i,j) = psi_i_store.*psi_j_store ;                  % Delta for store
           B(i,j) = (L-0)*trapz(y_L,psi_i_double.*psi_j_double) ;  % 
        end
    end
    % setup i and j torsion for wing and store, 2 torsion modes
    for i=1:M
        phi_i = (y_L).^i ;                                          % ith torsional function for wing
        phi_i_store = (Ls/L).^i ;                                   % ith torsional function for store
        phi_i_single = ((i)/(L)).*((y_L).^(i-1)) ;                   % ith torsional function first derivative (phi')
        for j = 1:M
            phi_j = (y_L).^j ;                                      % jth torsional function for wing 
            phi_j_store = (Ls/L).^j ;                               % jth torsional function for store
            phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;               % jth torsional function first derivative (phi')
            % D, D_store and T matrices
            D(i,j) = (L-0)*trapz(y_L, phi_i.*phi_j) ;                % D matrix for wing
            D_s(i,j) = phi_i_store.*phi_j_store ;                    % D matrix for store
            T(i,j) = (L-0)*trapz(y_L, phi_i_single.*phi_j_single) ;       % T matrix
        end
    end
    % Define functions for i bending and j torsion to determine C matrix
    for i=1:N
        psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
        psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
        for j = 1:M
            phi_j = (y_L).^j ;                                      % jth torsional function for wing
            phi_j_store = (Ls/L).^j ;                               % jth torsional function for store 
            phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;              % jth torsional function first derivative (phi')
            % C and C_store matrices
            C(i,j) = (L-0)*trapz(y_L, psi_i.*phi_j) ;
            C_s(i,j) = psi_i_store.*phi_j_store ;
        end
    end

% Expressions for wing and store mass matrix
Mwing = [   mbar*Del,   -mbar*x_b*C   ;
         -mbar*x_b*(C'),  Iwing*D   ] ;

Mstore = [    ms*Dels,     ms*xs*b*C_s   ; 
           ms*xs*b*(C_s'),   Is*D_s   ]   ;

%Determince expression for the stiffness matrix, unaffected by the addition of a store
K = [    EI*B,      zeros(N,M)  ;
     (zeros(N,M))',    GJ*T   ] ;

% Total mass matrix achieved by combining the mass matrices of store and wing
    Mt = Mwing + Mstore ;
    
for ii=1:150
    k=ii*0.01;
    % Calculate Theodorsen function
    C_theo=besselk(1,(j*k))./(besselk(0,(j*k))+besselk(1,j*k)) ;
    %calculate the A_mat, B_Mat, C_mat matrices
    A_mat = 2*pi*b*(k^2)*[Del a*b*C a*b*(C') (b^2)*((a^2)+(1/8))*D] ;
    B_mat = -2*pi*k*j*[2*C_theo]
    A_hat=A_mat+B_mat+C_mat;
    mu=eig(...);
    % Calculate omega from real part of mu
    omega(:,ii)=sqrt(1./real(mu));
    % Calculate damping from imaginary part of mu
    g(:,ii)=((omega(:,ii).^2).*imag(mu))/j;
    % calculate corresponding speed
    U(:,ii)=(omega(:,ii).*b)/k;
    % flutter offurs at zero damping
end

% check for solutions crossing the imaginary axis 
int=find(-0.5*imag(g(2,:))>0); 
int=int(end);
int=polyfit([U(2,int) U(2,(int+1))],[(-0.5*imag(g(2,int))) (-0.5*imag(g(2,(int+1))))],1);
Uf(int_pos)=abs(int(2))/abs(int(1));

end

% experimental dat afor mass "5"
Span_5=[0 0.23 0.375 0.54];
Uf_5=[ ...];

% data for mass "7e"
Span_7e=[0 0.23 0.44 0.6 0.73 0.83 1];
Uf_7e=[...];

% plots
...

% Repeat procedure for minimum flutter speed
% Determine the store position for minimum flutter speed
int_pos=find(Uf==min(Uf));
for ii=1:150
    k=ii*0.01;
    Ls=(int_pos-1).*(L/100);
    % Calculate Theodorsen function
    C_theo = bessel(1, (j*k))./(bessel(0, (j*k)) + bessel(1, j*k)) ;
    % Calculate Deta, Delta, B,C,D and T matrices
   
    y_L = 0:0.01:1 ;
    N = 3 ;
    M = 2 ;
    % calculate Delta, Deltas, B, C, D and T matrices setup functions for bending, 3 bending and 2 torsion modes 
   for i=1:N
        psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
        psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
        for j = 1:N
           psi_j = (y_L).^(j+1) ;                                   % jth bedning function of wing
           psi_j_store = (Ls/L).^(j+1) ;                            % jth bending mode of store
           psi_j_double = ((j*(j+1))/L.^2).*((y_L).^(j-1)) ;        % jth wing bending function second derivative (psi'') 
           % Del, Dels and B matrices
           Del(i,j) = (L-0)*trapz(y_L,psi_i.*psi_j) ;              % Delta for wing
           Dels(i,j) = psi_i_store.*psi_j_store ;                  % Delta for store
           B(i,j) = (L-0)*trapz(y_L,psi_i_double.*psi_j_double) ;  % 
        end
    end
    % setup i and j torsion for wing and store, 2 torsion modes
    for i=1:M
        phi_i = (y_L).^i ;                                          % ith torsional function for wing
        phi_i_store = (Ls/L).^i ;                                   % ith torsional function for store
        phi_i_single = ((i)/(L)).*((y_L).^(i-1)) ;                   % ith torsional function first derivative (phi')
        for j = 1:M
            phi_j = (y_L).^j ;                                      % jth torsional function for wing 
            phi_j_store = (Ls/L).^j ;                               % jth torsional function for store
            phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;               % jth torsional function first derivative (phi')
            % D, D_store and T matrices
            D(i,j) = (L-0)*trapz(y_L, phi_i.*phi_j) ;                % D matrix for wing
            D_s(i,j) = phi_i_store.*phi_j_store ;                    % D matrix for store
            T(i,j) = (L-0)*trapz(y_L, phi_i_single.*phi_j_single) ;       % T matrix
        end
    end
    % Define functions for i bending and j torsion to determine C matrix
    for i=1:N
        psi_i = (y_L).^(i+1) ;                                      % ith bending function for wing
        psi_i_store = (Ls/L).^(i+1) ;                               % ith bending funtion for store
        psi_i_double = ((i*(i+1))/L.^2).*((y_L).^(i-1)) ;           % ith bending function second derivative (psi'')
        for j = 1:M
            phi_j = (y_L).^j ;                                      % jth torsional function for wing
            phi_j_store = (Ls/L).^j ;                               % jth torsional function for store 
            phi_j_single = ((j)/(L)).*((y_L).^(j-1)) ;              % jth torsional function first derivative (phi')
            % C and C_store matrices
            C(i,j) = (L-0)*trapz(y_L, psi_i.*phi_j) ;
            C_s(i,j) = psi_i_store.*phi_j_store ;
        end
    end

% Expressions for wing and store mass matrix
Mwing = [   mbar*Del,   -mbar*x_b*C   ;
         -mbar*x_b*(C'),  Iwing*D   ] ;

Mstore = [    ms*Dels,     ms*xs*b*C_s   ; 
           ms*xs*b*(C_s'),   Is*D_s   ]   ;

%Determince expression for the stiffness matrix, unaffected by the addition of a store
K = [    EI*B,      zeros(N,M)  ;
     (zeros(N,M))',    GJ*T   ] ;

 
 % Determine expressions for wind and store mass matrix
    ...
%Determince expression for the stiffness matrix, unaffected by the
% addition of a store
    ...
% Total mass matrix achieved by combining the mass matrices of store and
% wing
    Mt=Mwing+Mstore;          
    
   ...
       
    A_hat=A_mat+B_mat+C_mat;
    mu= ...;
    % Calculate omega from real part of mu
    omega(:,ii)=sqrt(1./real(mu));
    % Calculate damping from imaginary part of mu
    g(:,ii)=((omega(:,ii).^2).*imag(mu))/j;
    % calculate corresponding speed
    U(:,ii)=(omega(:,ii).*b)/k; 
end

% plots 
...