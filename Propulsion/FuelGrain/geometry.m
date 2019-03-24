%% geometry
%{
    The purpose of this code is to analyze the inital and final diameter

    source - Stanford_AA284a_Lecture10
%}

clc;
clear;

Dp0 = 2.75:.05:4;
k = 10;

%paraffin N20
n = .5;
a = .155;

t = linspace(0,6.67,k);
dt = t(2)-t(1);
mox_dot = linspace(3.7,1.64,k+1);
Dp = reshape(repmat(0*t,1,length(Dp0)),[k,length(Dp0)])';
Dp = [Dp zeros(length(Dp0),1)];

figure(1);
hold on;


%Numeric Integration for circular grain port
for j = 1:length(Dp0)
    Dp(j,1) = Dp0(j);
    dDp_dt = (2 ^ (2*n + 1) * a / pi^n) * ( (mox_dot(1)^n) / (Dp(j,1)^(2*n)));
    for i = 1:length(t)
        Dp(j,i+1) = Dp(j,i) + dt*dDp_dt;
        dDp_dt = (2 ^ (2*n + 1) * a / pi^n) * ( (mox_dot(1+i)^n) / (Dp(j,1+i)^(2*n)));      
    end
    plot(t,Dp(j,2:length(Dp(1,:))));
end
%plot(t,Dp(1,2:length(Dp(1,:))));

L = 4*(17.78) / (pi*.036*(7.5)* (Dp(length(Dp0),11)^2 - Dp(length(Dp0),1)^2))
Dp(length(Dp0),1)
Dp(length(Dp0),11)



