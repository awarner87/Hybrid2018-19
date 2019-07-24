clear; clc; close all;
tic
% B.1
% Initial Conditions
%==========================================================================

% Tank Properties
%--------------------------------------------------------------------------
M_o   = 8.8843;                    % Initial tank fluid mass, kg (17.77 lbs)
V     = 0.0104386;                   % Tank volume, m3
T1_o  = 288.17;                  % Initial tank fluid temperature, K (70F)
%P_o  = 4.025;                    % Initial tank fluid pressure, Mpa
%X_o  = 0.01;                     % Initial tank fluid quality,
R     = 0.1889241;                % Gas constant, kJ/(kg*K)

% Injector Properties
%---------------------------------------------------------------------
d_inj = 0.125;       % Injector diameter, in (1.5-2mm x 13-15 holes), (0.21292625641732 in - 0.30495931858268 in)
Ac    = 16*(pi/4)*(d_inj*0.0254)^2;   % Injector cross sectional area; m3
Cd    = 0.65;                    % Injector discharge coefficient (0.65 - 0.75)
Pamb  = 2.75;                    % From Tank
Inj_Switch = 1;                  % Isentropic = 1, Adiabatic = 2
%Valve Specifications


% Time Iteration
%--------------------------------------------------------------------------
tstop = 300;                       % Stop time, s
dt    = 0.01;                       % Step time, s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize state
%--------------------------------------------------------------------------
rho1_o   = M_o/V;
Props1_o = N2Oprops(T1_o,rho1_o);
P1_o     = Props1_o.P;
X1_o     = Props1_o.X;
h1_o     = Props1_o.h;
H1_o     = M_o*Props1_o.h;
st_o     = Props1_o.state;

t    = 0;       % Column 1: Time (sec)
M    = M_o;     % Column 2: Tank fluid mass, M (kg)
rho1 = rho1_o;  % Column 3: Tank fluid density, rho1 (kg/m^3)
T1   = T1_o;    % Column 4: Tank temperature, T1 (K)
P1   = P1_o;    % Column 5: Tank pressure, P1 (MPa)
X1   = X1_o;    % Column 6: Tank quality, X1 ()
h1   = h1_o;    % Column 7: Tank specific enthalpy, h1 (kJ/kg)
H1   = H1_o;    % Column 8: Tank total enthalpy, H1 (J) ??? Units
mdot = 0;       % Column 9: Tank mass flow rate, mdot (kg/s)
P2   = Pamb;    % Column 10: Injector outlet pressure, P2 (Pa)
st   = st_o;    % Column 11: Fluid state, -1=-Input, 0=Liq, 1=Sat, 2=Gas
rhoL1 = Props1_o.rho_l; %Column 12: Liuid Density (kg/m^3)

State  = [t, M, rho1, T1, P1, X1, h1, H1, mdot, P2, st, rhoL1];
State2 = [0, 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 0;
while t < tstop
    i = i + 1;
    
    % Exit the loop when the tank is out of fluid
    %----------------------------------------------------------------------
    if M < 0; break; end
    
    % Initial Tank Fluid Properties
    %======================================================================
    %Note: this guess MUST yield a quality less than 1 and >0
    %guess=[300, 800]; %[T,rho]
%
    %Set up function (need to match pressure and quality)
    %pFunc = @(v) [getfield(N2Oprops(v(1),v(2)),'P')-P1; getfield(N2Oprops(v(1),v(2)),'X')-X1];
%{
    %Set up function (need to match pressure and quality)
    pFunc = @(v) {getfield(CO@PropsNIST(v(1),v(2)),'P')-P1 ...
                  getfield(CO2PropsNIST(v(1),v(2)),'X')-X1];
%}
    % Solve for [T,rho] of saturated but pure liquid at P1
    % lsqnonlin tries to find a T and rho that makes pFunc = 0
    Props1     = N2Oprops(T1, rho1);
%    Props1 = CO2PropsNIST(T1,rho1);
    Pv1     = Props1.P;                % Fluid Vapor Pressure, Mpa
    rhoL1  = Props1.rho_l;            % Fluid liquid density, kg/m3
    h1     = Props1.h;                % Fluid specific enthalpy, kJ/kg
    H1     = M*h1;                    % Fluid total enthalpy, kJ
    s1     = Props1.s;                %Fluid entropy, kJ/kg*K
    
    %Propagate properties across the injector
    %======================================================================
    %Assume an isentropic (s1=s2) or adiabatic (h1=h2) injector, solve for properties at P2
    if Inj_Switch == 1
        %Set up function to match pressure and entropy
 %
        pFunc = @(v) [getfield(N2Oprops(v(1),v(2)),'P')-P2 ...
                      getfield(N2Oprops(v(1),v(2)),'s')-s1];

%{
        pFunc = @(v) [getfield(CO2PropsNIST(v(1),v(2)),'P')-P2 ...
                      getfield(CO2PropsNIST(v(1),v(2)),'s')-s1];
%}
    elseif Inj_Switch == 2
        %Set up function to match pressure and enthalpy
        pFunc = @(v) [getfield(N2Oprops(v(1),v(2)),'P')-P2 ...
                      getfield(N2Oprops(v(1),v(2)),'h')-h1];
    end
    
    % Solve for T2 & rho2 downstream of the injector
    %----------------------------------------------------------------------
    guess=[300, 800]; %[T, rho]
    v2 = ...
     lsqnonlin(pFunc,guess,[0, 0],[inf, inf],optimset('Display','off','TolFun',1e-14));
    T2 = v2(1); rho2 = v2(2);
     
    % Solve for the remainng properties downstream of the injector
    %----------------------------------------------------------------------
    Props2=N2Oprops(T2,rho2);          % Downstream fluid properties
%     Props2=CO2PropsNIST(T2,rho2);

    Pv2   = Props2.P;                  % Downstream fluid vapor pressure
    h2    = Props2.h;                  % Downstream fluid enthalpy
    
    % Calculate the Mass Flow
    %======================================================================
    % Non-Equilibrium Parameter
    %k = sqrt(abs(P1-P2)/abs(Pv2-P1)); %Whitmores Equation
    k = sqrt((P1-P2)/(Pv1-P2));       % Spencer and Stanfords Equation
    
    % Weighting Coefficeint
    W = (1/(k+1));
    
    % Incompressible fluid mass flow rate
    mdot_inc = Ac*sqrt(2*rhoL1*(P1-P2)*1e6);
    if imag(mdot_inc) ~= 0; break; end
    
    % Homogeneous Equilibrium mass flow rat
    mdot_HEM = rho2*Ac*sqrt(2*(h1-h2));
    if imag(mdot_HEM) ~= 0; break; end

    %Weighted Non-Homogenous Equilibrium (modified Standford) mass flow
    %rate
    mdot = Cd*((1-W) * mdot_inc + W * mdot_HEM); % Shannon's Theory
    
    % Update the upstream fluid properties for the nest step
    %======================================================================
    M    = M - mdot*dt;              % Update tank fluid mass
    if M < 0; break; end
    Hdot = h1*mdot;                   % Enthalpy flow rate
    H1   = H1 - Hdot*dt;              % Update tank total enthalpy
    
    % Calculate the new tank enthlapy and density
    %----------------------------------------------------------------------
    rho1 = M/V;                      % Update tank specific density
    h1   = H1/M;                      % Update tank specific enthalpy
    
    % Calculate the new tank temperature
    %----------------------------------------------------------------------
    % Create a function for lsqnonlin to solve T(rho,h)
    pFunc = @(T_Unknown) getfield(N2Oprops(T_Unknown,rho1),'h')-h1;
%     pFunc = @(T_Unknown) getfield(CO2PropsNIST(T_Unknown,rho1),'h')-h1;
    %Since T_Unknown is not pre define in pFunc, lsqnonlin will fin a T for rho_Known and h_Known
    T1=lsqnonlin(pFunc,300,0, inf,optimset('Display','off','TolFun',1e-14));
%      [T1] = CO2_rho_h_2T(rho1,h1,300);

    %Calculate the new tank pressure and quality
    %----------------------------------------------------------------------
    Props1 = N2Oprops(T1,rho1);
%     Props 1 = CO2PropsNIST(T1,rho1);
    P1 = Props1.P;
    X1 = Props1.X;
    st1 = Props1.state;
    
    % Update the state
    %----------------------------------------------------------------------
    t = t + dt;
    
    State = [State; t, M, rho1, T1, P1, X1, h1, H1, mdot, Pamb, st1, rhoL1];
    State2 = [State2; mdot_inc, mdot_HEM];
end
%
% Open cold flow data for comparison
%==========================================================================
%%%%%%CF = open('Cold Flow Data\FillTest1_29_largeNoTop.mat');
%CF = open('Cold Flow Data\FillTest1_29_smallNoTop.mat');
%CF = open('Cold Flow Data\FillTest1_29_large.mat');
%CF = open('Cold Flow Data\FillTest1_29_small.mat');

% % % % Find portion of cold flow data where the run valve is open
% % % %--------------------------------------------------------------------------
% % % iter = 0;
% % % for loop = 1:length(CF.Valvedat)
% % %     if CF.Valvedat(loop,5) == 1
% % %         iter = iter + 1;
% % %         CF_t(iter)           = CF.t(loop);
% % %         CF_Weight(iter)      = CF.Weight(loop);
% % %         CF_PTank(iter)       = CF.PRunTank(loop);
% % %         CF_Mdot(iter)        = CF.mdotVentIn(loop);
% % %         CF_DeltaWeight(iter) = CF.DeltaWeightKg(loop);
% % %         CF_T_Tank(iter)      = CF.Tank_temp(loop);
% % %         CF_dP_vent(iter)     = CF.dpVenturi(loop);
% % %         CF_rho1_H(iter)      = CF.rho1_H(loop);
% % %     end
% % % end
% % % CF_t = CF_t - CF_t(1);

% Calculate the venturi mass flow rate the helmholts density
%--------------------------------------------------------------------------
d1        = .0254762;
d2        = .0078486;
A1        = pi*d1^2/4;
A2        = pi*d2^2/4;
C         = 0.980243097172394;
%%%CF_Mdot_H = A2.*C.*CF_rho1_H.*sqrt(2.*(CF_dP_Vent)./(CF_rho1_H.*(1-(A2/A1)^2)));

% % % % Calculate the mass flow rate by differentiating the tank mass
% % % %--------------------------------------------------------------------------
% % % CF_Mdot_diff = - diff(CF_Weight)./diff(CF_t);
% % % CF_dt        = CF_t(2) - CF_t(1);
% % % %CF_Mdot_wave = derivative_cwt(CF_Weight,'gaus1',16,CF_dt,1);
% % % CF_Mdot_wave = derivative_dwt(CF_Weight,'sp1',4,CF_dt,1);

% % % % Open cold flow data for comparison
% % % %==========================================================================
% % % Whitmore = open('Whitmores Output\FillTest1_29_largeNoTop.mat');
% % % W_t    = Whitmore.t;
% % % W_M    = Whitmore.M;
% % % W_rho1 = Whitmore.rho1;
% % % W_T1   = Whitmore.T1;
% % % W_P1   = Whitmore.P1;
% % % W_X1   = Whitmore.X1;
% % % W_H1   = Whitmore.H1;
% % % W_mdot = Whitmore.mdot;
% % % W_2    = Whitmore.P2;

% Up-pack the state vector for plotting
%==========================================================================
t    = State(:,1); % Column 1: Time (sec)
M    = State(:,2); % Column 2: Tank fluid mass, M (kg)
rho1 = State(:,3); % Column 3: Tank fluid density, rho1 (kg/m^3)
T1   = State(:,4); % Column 4: Tank temperature, T1 (K)
P1   = State(:,5); % Column 5: Tank pressure, P1 (Pa)
X1   = State(:,6); % Column 6: Tank quality, X1 ()
h1   = State(:,7); % Column 7: Tank specific enthalpy, h1 (kJ/kg)
H1   = State(:,8); % Column 8: Tank total enthalpy, H1 (J) ??? Units
mdot = State(:,9); % Column 9 : Tank mass flow rate, mdot (kg/s)
P2   = State(:,10);% Column 10: Injector outlet pressure, P2 (Pa)
st   = State(:,11);% Column 11: Fluid state, -1=-Input, 0=Liq, 1=Sat, 2=Gas
rhoL1 = State(:, 12); %Volumn 12: Liquid Density (kg/m^3)

mdot_inc = State2(:,1); % Mdot Incompressible
mdot_HEM = State2(:,2); % Mdot HEM

%%Plot
%==========================================================================
% Column 2 - Tank fluid mass
%--------------------------------------------------------------------------
figure(2)
plot(t,M); hold on; grid on;
%plot(W_t,W_M,'--g');
%%%plot(CF_t,CF_Weight,'--r');
%plot(CF_t,CF_DeltaWeight,'g');
title('Tank Mass');
%xlim([0 tstop]);
xlabel('Time, s'); ylabel('Tank Fluid Mass, kg');
%legend('Model Prdicted', 'Whitmore', 'ColdFlow');
%legend('Model Predicted', 'Cold Flow');
%%%saveas(2,'Tank_Mass.Png')

% Column 3 - Tank fluid specific density
%--------------------------------------------------------------------------
figure(3);
plot(t,rho1); hold on; grid on;
%%%plot(W_t,W_rho1,'--g');
%%%plot(CF_t,CF_Weight/V,'r');
title('rho1');
xlabel('Time, s'); ylabel('Density, kg/m3');
%legend('Model Predicted', 'Whitmore', 'Cold Flow');
%%%saveas(3, 'Tank_Density.Png')

% Column 4 - Tank temperature
%--------------------------------------------------------------------------
figure(4);
plot(t,T1); hold on; grid on;
%plot(W_t,W_T1,'--g');
%%%plot(CF_t,CF_T_Tank,'--r');
title('Temperature');
%xlim([0 tstop]);
xlabel('Time, s'); ylabel('Temperature, K');
%legen('Model Predicted', 'Whitmore', 'Cold Flow');
%legend('Model Predicted','Cold Flow');
%%%saveas(4,'Tank_Temperature.Png')

% Column 5 - Tank pressure
%--------------------------------------------------------------------------
figure(5)
plot(t,P1); hold on; grid on;
%plot(W_t,W_P1/1e3,'--g');
%%%plot(CF_t,CF_PTank/1e6,'--r');
title('P1');
%xlim([0 tstop]);
xlabel('Time, s'); ylabel('Pressure, MPa');
%legend('Model Predicted', 'Cold Flow');
%legend('Model Predicted', 'Whitmore', 'Cold Flow');
%%%saveas(5,'Tank_Pressure.Png')

% Column 6 - Tank quality
%--------------------------------------------------------------------------
figure(6)
plot(t,X1); hold on; grid on;
%%%plot(W_t,W_X1,'--g');
ylim([0 1]);
title('X1');
xlabel('Time,s'); ylabel('Quality');
%legend('Model Predicted','Whitmore');
%%%saveas(6,'Tank_Quality.Png')

% Column 7 - Tank specific enthalpy
%--------------------------------------------------------------------------
figure(7)
plot(t,h1); hold on; grid on;
title('h1');
xlabel('Time, s'); ylabel('Specific Enthalpy, kJ/kg');
%legend('Model Predicted');
%%%saveas(7,'Tank_Specific_Enthalpy.Png')

% Column 8 - Tank total enthlapy
%--------------------------------------------------------------------------
figure(8)
plot(t,H1); hold on; grid on;
%%%plot(W_t,W_H1,'--g');
title('H1'); xlabel('time,s'); ylabel('Enthalpy, kJ');
%legend('Model Predicted', 'Whitmore');
%%%saveas(8,'Tank_Total_Enthalpy.Png')

% Column 9 - Mass flow rate 
%--------------------------------------------------------------------------
figure(9)
plot(t,mdot,'b'); hold on; grid on;
% % % plot(CF_t,CF_Mdot,'--r')
% % % plot(CF_t,CF_Mdot_H,'-.b')
title("Model Mass Flow Rate"); %%%THIS LINE WAS ADDED
xlabel('Time,s'); ylabel('Mass Flow Rate, kg/s');
%xlim([0 tstop]); %ylim([0 3]);
%legend('Model Predicted','Venturi-Curve Fit Denstiy','Venturi-Helmholtz Density');
%%%saveas(9,'Tank_Mass_Flow_Rate_Venturi.Png')

figure(90);
plot(t,mdot,'b'); hold on; grid on;
%plot(W_t,W_mdot,'--g');
% % % plot(CF_t,CF_Mdot,'--r')
% % % plot(CF_t,-CF_Mdot_wave,'-.b')
plot(t,mdot_inc,':b');
plot(t,mdot_HEM,'--bo');
title("Mass Flow Rate Breakdown") %%%THIS LINE WAS ADDED
xlabel('Time, s'); ylabel('Mass Flow Rate, kg/s');
%xlim([0 tstop]); ylim([0 1.6]);
legend('Model Predicted','Incompressible Model','HEM Model');
%%%saveas(90,'Tank_Mass_Flow_Rate.Png')

% % % figure(900);
% % % plot(CF_t(2:end),CF_Mdot_diff); hold on; grid on;
% % % plot(CF_t,-CF_Mdot_wave,'--r')
% % % xlabel('Time, s'); ylabel('Mass Flow Rate, kg/s');
% % % xlim([0 tstop]); ylim([-2 4]);
% % % legend('Cold Flow-dM/dt-diff','Cold Flow-dM/dt-wavelet');
% % % saveas(900,'Tank_Mass_Flow_Rate_Diff.Png')

% Column 11 - Tank fluid state
%--------------------------------------------------------------------------
figure(11);
plot(t,st); grid on;
title('Fluid State');
xlabel('Time, s'); ylabel('Fluid State');
ylim([-2 3]);
%legend('Model Predicted');
%%%saveas(11,'Tank_Fluid_State.Png')
%
toc

mdot_max = max(mdot)*2.20462;
mdot_avg = mean([mdot(2:end-1); 0])*2.20462;
fprintf("Maximum Mass Flow Rate: %.2f lb/s\n", mdot_max)
fprintf("Average Mass Flow Rate: %.2f lb/s\n", mdot_avg)
fprintf("Time to Tank Emptry: %.2f s\n", t(end))

%% Pressure drop calculation (Needs Improvement)
alpha = 1;
pipe_D = .402; %in
A_pipe = (pipe_D * .0254)^2 * pi; %m^2
rho_pipe = rho1;
V_A = mdot ./ (rho1*A_pipe); %m/s
g = 7*9.8; %7g, m/s
P_A = P1;

%A to B (Tank to Main Valve Entrance)
L_AB = 2; %in
delZ_AB = 2; %in
f = .02;
h_LT_AB = f * L_AB/pipe_D * V_A.^2/(2*g); %m
P_B = P_A + (rho_pipe .* g .* (delZ_AB .* .0254 - h_LT_AB))/1000000; %MPa

%B to C (Main Valve Entrance to Main Valve Exit)
cv = 12;
q = mdot./rho_pipe*15850.372483753; %GPM %.5 * rho_pipe .* V_A.^2; %Pa
SG = rho_pipe/997;
delP_BC = q.^2./cv^2 .* SG; %psi
P_C = P_B - delP_BC/145.038; %MPa

%C to D (Main Valve Exit to Injector Plate)
L_CD = 2; %in
delZ_CD = 2; %in
f = .02;
h_LT_CD = f * L_CD/pipe_D * V_A.^2/(2*g); %m
P_D = P_C + (rho_pipe .* g .* (delZ_CD .* .0254 - h_LT_CD))/1000000; %MPa

%Calculate and Plot Pressure Difference
delP_AD = (P_A - P_D)*145.038; %psi
delP_D2 = (P_D - P2)*145.038; %psi
figure;
hold on;
plot(t, delP_AD);
plot(t, delP_D2);
title("Pressure Drop Comparison");
ylabel("Pressure Drop (psi)");
xlabel("Time (s)");
legend("Tank to Injector Plate", "Across Injector Plate");
axis([0, max(t), min([zeros(size(delP_AD)); delP_AD; delP_D2]), max([delP_AD; delP_D2])]);
grid on

%Plot Pipe Velocity

V_pipe_comp = mdot ./ (rho1 .* A_pipe);
figure;
hold on
plot(t, V_A);
plot(t, V_pipe_comp);
title("Velocity from Tank to Injector");
xlabel("Time (s)");
ylabel("Velocity (m/s)");

grid on
