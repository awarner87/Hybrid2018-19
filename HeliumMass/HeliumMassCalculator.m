%OBTAINED FROM https://books.google.com/books?id=TKdIbLX51NQC&pg=PA141&lpg=PA141&dq=calculating+pressurant+mass+in+a+rocket&source=bl&ots=sma25Mguu3&sig=ACfU3U0wyUI0mb7tuySXm1mowbsXbDwRJw&hl=en&sa=X&ved=2ahUKEwi13qfIxajkAhUEGKwKHcTTAQIQ6AEwAXoECAgQAQ#v=onepage&q=calculating%20pressurant%20mass%20in%20a%20rocket&f=false
%%ALL DONE ASSUMING PRESSURE OF 1000 psi%%
Data_N2O = xlsread('N2O50.100.xlsx');
Data_He = xlsread('Helium50.100.xlsx');
M_N2O = 17.78;              %lbm
T_N2O = Data_N2O(:,1);      %F
T_He = Data_He(:,1);        %F
Rho_N2O = Data_N2O(:,3);    %lbm/ft^3
Rho_He = Data_He(:,3);      %lbm/ft^3
P_N2O = 1000;               %psia
P_He = 1000;                %psia
V_total = 782/(12^3)*ones([numel(T_N2O),1]);   %ft^3
V_N2O = M_N2O./Rho_N2O;     %ft^3
V_He = V_total - V_N2O;     %ft^3
M_He = V_He.* Rho_He; %lbm
initialM_He = V_total.*Rho_He;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Helium Requirements');
subplot(3,1,1);
plot(T_N2O,initialM_He,'b+');
hold on;
grid on;
grid minor;
xlabel('Temperature of Helium (F)');
ylabel('Mass of Helium (lbm)');
title('Mass of Helium Required to Initially Fill Ox Tank To 1000 Psi');
datacursormode;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2);
plot(T_N2O, M_He,'r+');
hold on;
grid on;
grid minor;
xlabel('Temperature of Fluids (F)');
ylabel('Mass of Helium (lbm)');
title('Mass of Helium Required to Achieve 1000 psi in Ox Tank');
datacursormode;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,3);
plot(T_N2O, V_He,'k+');
hold on;
grid on;
grid minor;
xlabel('Temperature of Fluids (F)');
ylabel('Volume of Helium (ft^3)');
title('Volume of Helium Required to Achieve 1000 psi in Ox Tank');
datacursormode;
hold off;