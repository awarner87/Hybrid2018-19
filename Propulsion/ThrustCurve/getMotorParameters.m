function [oxFlux,regRate,grainMdot,OtoF,chambPress,motorThrust] = getMotorParameters(oxMdot,grainDport,grainLength,nozzleArea)
%Author: Austin Keck, keck2@purdue.edu
%Last Updated: 10/29/18
%Dependent Scripts: None

%DESCRIPTION%
%Based on Stanford lecture notes: http://www.spg-corp.com/docs/Stanford_AA284a_Lecture10.pdf
%Calculates instantaneous motor parameters for getThrustCurve.

%ASSUMPTIONS/NOTES%
%Space-averaged port diameter
%Circular port
%Complete combustion
%Chamber pressure, thrust, impulse calculations incomplete

%INPUTS% 
%oxMdot: oxidizer mass flow rate (lb/sec)
%grainDport: initial solid grain port diameter (in)
%grainLength: solid grain length (in)
%nozzleArea: nozzle throat area (in^2)

%OUTPUTS%
%oxFlux: instantaneous oxidizer mass flux (lb/sec*in^2)
%regRate: instantaneous regression rate (in/sec)
%OtoF: instantaneous O/F ratio (-)
%chambPress: instantaneous chamber pressure (psi)
%motorThrust: instantaneous motor thrust (lbf*ft/sec^2)

%DEMO-INPUTS%
%{
Demo 1: 
[oxFlux,regRate,grainMdot,OtoF,chambPress,motorThrust] = getMotorParameters(2.67,16,1)
%}

%INITIALIZATION%
regCoefA = 0.161807; %paraffin regression coefficient
regExpN = 0.5; %paraffin regression exponent
rohParaf = 0.032514563; %paraffin density (lb/in^3)
cStarTheo = 205; %theoretical characteristic velocity
etaC = 0.85; %Efficiency of combustion from Stanford example (TBD)
etaN = 0.97; %Efficiency of nozzle from Stanford example 
coefThrust = 0.77; %(TBD)

oxFlux = 4 * oxMdot / (pi * grainDport .^ 2);
regRate = regCoefA * oxFlux .^ regExpN;
grainMdot = rohParaf * pi * grainDport * grainLength * regRate;
OtoF = oxMdot / grainMdot;

chambPress = (oxMdot + grainMdot) * cStarTheo * etaC / nozzleArea;
motorThrust = (oxMdot + grainMdot) * cStarTheo * etaC * coefThrust * etaN;