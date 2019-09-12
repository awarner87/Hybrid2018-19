%function [W_He] = HeliumMassCalculator(p1,p2,T1);
%OBTAINED FROM https://books.google.com/books?id=TKdIbLX51NQC&pg=PA141&lpg=PA141&dq=calculating+pressurant+mass+in+a+rocket&source=bl&ots=sma25Mguu3&sig=ACfU3U0wyUI0mb7tuySXm1mowbsXbDwRJw&hl=en&sa=X&ved=2ahUKEwi13qfIxajkAhUEGKwKHcTTAQIQ6AEwAXoECAgQAQ#v=onepage&q=calculating%20pressurant%20mass%20in%20a%20rocket&f=false
%%ALL DONE ASSUMING PRESSURE OF 1000 psi%%
M_N2O = 17.78; %lbm
T_N2O = (530:1:553)';
Rho_N2O = [50.37 50.094 49.811 49.521 49.223 48.916 48.601 48.275 47.938 47.589 47.226 46.848 46.452 46.036 45.597 45.131 44.632 44.093 43.505 42.852 42.110 41.236 40.125 38.581]';    %lbm/ft^3
P_N2O = 1000;   %psi
V_total = 782/12^3*ones([numel(T_N2O),1]);   %in^3
V_N2O = M_N2O./Rho_N2O;  %ft^3
V_He = V_total - V_N2O;  %ft^3
Rho_He = [0.68170 0.68046 0.67923 0.678 0.67677 0.67555 0.67434 0.67312 0.67192 0.67071 0.66951 0.66832 0.66713 0.66594 0.66476 0.66358 0.66241 0.66124 0.66007 0.65891 0.65775 0.6566 0.65545 0.65431]';    %lbm/ft^3

%R_N2O = 35.112;     %ft-lbf/lb-R
%Z = P_N2O./Rho_N2O/R_N2O./T_N2O;

%NEGLECTING ALL HEAT TRANSFER% ]
M_He = V_He.* Rho_He; %lbm
figure;
plot(T_N2O,M_He,'r*');
hold on;
grid on;
grid minor;
xlabel('Temperature of Fluids (R)');
ylabel('Mass of Helium (lbm)');
title('Mass of Helium Required to Achieve 1000 psi in Ox Tank');
datacursormode;
hold off;

figure;
plot(T_N2O,V_He,'g*');
hold on;
grid on;
grid minor;
xlabel('Temperature of Fluids (R)');
ylabel('Volume of Helium (ft^3)');
title('Volume of Helium Required to Achieve 1000 psi in Ox Tank');
datacursormode;
hold off;
%end