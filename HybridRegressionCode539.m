% Madeleine Melcer
% AAE 339 HW 5 Question 1

clear all
close all

%% given values
w_a1 = 0.6; % web fraction 
w_a2 = 0.8;
w_b = 0.65;
rho = 1760; % kg/m^3
Cstar = 1568; % m/s
Tc = 3392; % K
gam = 1.16; 
a = 0.561/100; % m/s-MPa^n
n = 0.35;
Mwg = 22; % g/mol
R_u = 8.314; % J/mol-K
g = 9.81; % m/s^2
Pc_initial = 5; % MPa

% conversions
lb2kg = 1/2.205; % conversion from lbm to kg

%% part i

LD_ratio = @(w) 1/(2*w)*(3-(1-w)^2-2*(1-w));

LD_a1 = LD_ratio(w_a1);
LD_a2 = LD_ratio(w_a2);

fprintf('Part i: \nL/Do for web fraction of %.1f is %.1f \n',w_a1,LD_a1);
fprintf('L/Do for web fraction of %.1f is %.1f \n\n',w_a2,LD_a2);

%% part ii

LD_b = LD_ratio(w_b);

mp = 1*lb2kg; % kg, given that there is 1 lb of propellant in this motor
V = mp/rho; % kg/m^3
ro = (V/(2*LD_b*pi*(1-(1-w_b)^2)))^(1/3); % m
ri = (1-w_b)*ro; % m
web = ro-ri; % m
D = 2*ro; % m
L = LD_b*D; % m

fprintf('Part ii: \nL is %.3f cm, Do is %.3f cm, and web thickness is %.3f cm \n\n',...
    L*100,D*100,web*100);

%% part iii

A_b = @(w) 2*pi*(ri+w)*(L-2*w) + 2*pi*(ro^2-(ri+w)^2); % m^2
A_b_initial = A_b(0); % m^2

A_t = a*rho*A_b_initial*Cstar/Pc_initial^(1-n)/(10^6); % m^2
D_t = 2*sqrt(A_t/pi); % cm

fprintf('Part iii: \nThe throat diameter is %.3f cm \n\n',D_t*100);

%% part iv

dw = 0.001; % m

i = 1;
w(1) = 0; % m
Ab(1) = A_b_initial; % m^2
Pc(1) = 5; % MPa
r(1) = a*Pc(1)^n; % m/s
t(1) = 0; % s

while w(i) < web
    i = 1 + i;
    w(i) = w(i-1) + dw; % m
    Ab(i) = A_b(w(i)); % m^2
    Pc(i) = (a*rho*Ab(i)*Cstar/A_t/(10^6))^(1/(1-n)); % MPa
    r(i) = a*Pc(i)^n; % m/s
    t(i) = t(i-1) + dw/r(i); % s
end

R = R_u/Mwg*1000; % J/kg-K
rho_g = Pc_initial*(10^6)/R/Tc; % kg/m^3
mdot = Pc_initial*(10^6)*A_t/Cstar; % kg/s
tau = rho_g*(pi*ro^2*(L+2*web))/mdot; % s

i = 1;
t_blowdown(1) = t(end); % s
Pc_blowdown(1) = Pc(end); % MPa
while Pc_blowdown(i)>0.0001
    i = i+1;
    t_blowdown(i) = t_blowdown(i-1) + 0.001;
    Pc_blowdown(i) = Pc(end)*exp(-(t_blowdown(i)-t(end))/tau); % MPa
end

figure(1)
hold on
plot(t,Pc)
plot(t_blowdown,Pc_blowdown)
xlabel('Time (s)')
ylabel('Pc (MPa)')
title('Question 1: Pressure history of motor')
legend('normal operation','blowdown')
hold off
