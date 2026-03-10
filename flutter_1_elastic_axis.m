% Hard Flutter Solution for wing with mass-balanced section (i.e. elastic axis and center of mass coinside xa=0)

clear 
clc

ra = 0.5 ;          % radius of gyration
mu = 10 ;           % mass ratio
fr = 0.5 ;          % frequency ratio omega_h / omega_a
formatSpec = 'UFhat_F %f, wFwa %f, nu %f\n';

for jj = 1:71
E = 0.12 + (jj - 1)*.005 ;      % Elastic axis ec/b=0.5 or default E=0.25
UFhat = 0.5 ;                   % guess reduced flutter speed
wFwa = 0.1 ;                    % guess flutter frequency
nu = 2*wFwa/UFhat ;             % reduced flutter frequency
fprintf (formatSpec,UFhat,wFwa,nu)
    for ii = 1:20
    CF=besselk(1,(j*nu/2))./(besselk(0,(j*nu/2))+besselk(1,j*nu/2)) ;
    Muthetadot=2*pi*(-0.25*(0.5-E)+E*real(CF)*(0.5-E)+E*imag(CF)/nu) ;
    wFwa=sqrt((1-fr*fr*2*Muthetadot/pi/ra/ra)/(1-2*Muthetadot/pi/ra/ra)) ;
    nu=2*wFwa/UFhat ;
    UFhat=sqrt(0.5*mu*ra*ra*( (wFwa*wFwa*(wFwa*wFwa-1-fr*fr)+fr*fr)/(2*fr*fr*E-wFwa*wFwa*(2*E+4*Muthetadot/mu/pi))));  
    fprintf (formatSpec,UFhat,wFwa,nu)
    end
    fprintf (formatSpec,UFhat,wFwa,nu)
    UFHAT(jj)=UFhat;
    EE(jj)=E;
end

figure('Name', 'Flutter - Elastic Axis')
plot (EE, UFHAT, 'o', 'LineWidth', 2)
title('Effect of Elastic Axis Placement on Flutter Speed')
xlabel('Elastic axis location (ec/b)')
ylabel('Reduced flutter speed (U_F/(b\omega_a))');
grid on
