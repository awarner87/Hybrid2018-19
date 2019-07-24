function [props] = N2Oprops(T, rho)
%Returns structure with properties of N2O
% Helmholtz: Valid for T < 555 K and p < 50 MPa
% Saturation: Valid for 183.15 K < T < 309.15 K
% Inputs:
%   Temperature T (K)
%   Density rho (kg/m^3)
% Outputs:
%   Pressure P (MPa)
%   Specific Internal Energy u (kJ/kg)
%   Specific Enthalpy h (kJ/kg)
%   Specifc Entropy s (kJ/(kg*K))

%% CONSTANTS
Tc = 309.52; %Critical Temperature, K
rhoc = 10.27*44.0128; %Critical Density, kg/m^3
pc = 7.245; %Critical Pressure, MPa
R = 8.3144598/44.0128*1000; %Specific Gas Constant, kJ/(kg*K);

%% Nondimensionalize Inputs
t = Tc/T; %Inverse reduced temperature, dimensionless

%% Saturation Properties
bl = [1.72328, -0.83950, 0.51060, -0.10412];
bg = [-1.00900, -6.28792, 7.50332, -7.90463, 0.629427];
rhoL = rhoc * exp(bl(1)*(1-1/t)^(1/3) + bl(2)*(1-1/t)^(2/3) + bl(3)*(1-1/t) + bl(4)*(1-1/t)^4/3);
rhoG = rhoc * exp(bg(1)*(t-1)^(1/3) + bg(2)*(t-1)^(2/3) + bg(3)*(t-1) + bg(4)*(t-1)^(4/3) + bg(5)*(t-1)^(5/3));
bp = [-6.71893, 1.35966, -1.3779, -4.051];
pSat = pc*10^6*exp(t*(bp(1)*(1-1/t) + bp(2)*(1-1/t)^(3/2) + bp(3)*(1-1/t)^(5/2) + bp(4)*(1-1/t)^(5)));

%% CHECK STATE
if T > Tc
    state = 3; %Supercritical
    X = NaN;
    [p, u, h, s] = Helmholtz(T, rho);
else
    if rho >= rhoL
        state = 0; %Compressed Liquid
         X = NaN;
        [p, u, h, s] = Helmholtz(T, rho);
    elseif rho <= rhoG
        state = 2; %Superheated Vapor
        X = NaN;
        [p, u, h, s] = Helmholtz(T, rho);
    else
        state = 1; %Saturated Liquid Vapor Mix
        X = rhoG/rho * (rho - rhoL)/(rhoG - rhoL);
        [pl, ul, hl, sl] = Helmholtz(T, rhoL);
        [pg, ug, hg, sg] = Helmholtz(T, rhoG);
        p = pSat;
        u = ul + X*(ug - ul);%ug*X + ul*(1-X);
        h = hl + X*(hg - hl);%hg*X + hl*(1-X);
        s = sl + X*(sg - sl);%sg*X + sl*(1-X);
    end
end


%% CREATE OUTPUT STRUCTURE
props.P = p/10^6;
props.u = u/1000; %kJ/kg
props.h = h/1000; %kJ/kg
props.s = s/1000; %kJ/(kg*K)
props.state = state;
props.X = X;
props.rho_l = rhoL;
props.rho_v = rhoG;

end
function [p, u, h, s] = Helmholtz(T, rho)
%% CONSTANTS
Tc = 309.52; %Critical Temperature, K
rhoc = 10.27*44.0128; %Critical Density, kg/m^3
R = 8.3144598/44.0128*1000; %Specific Gas Constant, kJ/(kg*K);

%% Nondimensionalize Inputs
t = Tc/T; %Inverse reduced temperature, dimensionless
d = rho/rhoc; %Reduced density, dimensionless

%% IDEAL GAS PART
v = [2.1769, 1.6145, 0.48393]; %Coefficients for cp0, dimensionless
u = [879, 2372, 5447]; %Coefficients for cp0, K
c0 = 3.5; %Constant term for cp0
a = [-4.4262736272, 4.3120475243]; %Coefficients for ideal gas Helmholtz
a0 = a(1) + a(2)*t + log(d)+ (c0-1)*log(t) + sum(v .* log(1-exp(-u.*t/Tc))); %Dimensionless Helmholtz Energy, Ideal Gas Part

%% RESIDUAL PART
i = [1, 1, 1, 3, 7, 1, 2, 5, 1, 1, 4, 2]; %Reduced density exponent, same for all fluids
j = [0.25, 1.25, 1.5, 0.25, 0.875, 2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5]; %Inverse reduced temperature exponent, same for all fluids
l = [1, 1, 1, 2, 2, 2, 3]; %Reduced density exponential exponent, same for all fluids
k = [0.88045, -2.4235, 0.38237, 0.068917, 0.00020367, 0.13122, 0.46032, -0.0036985, -0.23263, -0.00042859, -0.042810, -0.023038]; %Coefficients for Equation of State, N2O Specific
ar = sum(k(1:5).*d.^i(1:5).*t.^j(1:5)) + sum(k(6:12).*d.^i(6:12).*t.^j(6:12).*exp(-d.^l)); %Dimensionless Helmholtz Energy, Residual Part

%% Derivative Calculation
a0_t = a(2) + (c0-1)/t + sum(u.*v.*exp(-u.*t/Tc)./(Tc.*(1-exp(-u.*t/Tc))));
ar_t = sum(k(1:5).*d.^i(1:5).*j(1:5).*t.^(j(1:5)-1)) + sum(k(6:12).*d.^i(6:12).*j(6:12).*t.^(j(6:12)-1).*exp(-d.^l));
ar_d = sum(k(1:5).*i(1:5).*d.^(i(1:5)-1).*t.^j(1:5)) + sum(k(6:12).*t.^j(6:12).*(i(6:12)-l.*d.^l).*d.^(i(6:12)-1).*exp(-d.^l));

%% FLUID PROPERTY CALCULATIONS
Z = 1 + d * ar_d; %Compressibility Factor, dimensionless
p = rho * R * T * Z; %Pressure, Pa
u = R*T*t*(a0_t+ar_t); %Specific Internal Energy, J/kg
h = u + R*T*(d*ar_d + 1); %Specific Enthalpy, J/kg
s = R*(t*(a0_t + ar_t) - a0 - ar); %Specifc Entropy, J/(kg*K)

end

