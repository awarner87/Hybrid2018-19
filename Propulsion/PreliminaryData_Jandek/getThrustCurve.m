function [totalImp,grainDport,meanOtoF] = getThrustCurve(timeHist,oxMdotHist,grainDport,grainLength,nozzleArea)
%Author: Austin Keck, keck2@purdue.edu
%Last Updated: 10/29/18
%Dependent Scripts: getMotorParameters

%DESCRIPTION%
%Based on Stanford lecture notes: http://www.spg-corp.com/docs/Stanford_AA284a_Lecture10.pdf
%Models transient solid fuel grain regression rate, O/F, and impulse given
%combustion chamber geometry and oxidizer flow history.

%ASSUMPTIONS/NOTES%
%Space-averaged port diameter
%Circular port
%Complete combustion
%Chamber pressure, thrust, impulse calculations incomplete

%INPUTS% 
%timeHist: time history data (sec)
%oxMdotHist: oxidizer mass flow rate (lb/sec)
%grainDport: initial solid grain port diameter (in)
%grainLength: solid grain length (in)
%nozzleArea: nozzle throat area (in^2)

%OUTPUTS%
%totalImp: total impulse (lbf*sec)
%grainDport: final solid grain port diameter (in)
%meanOtoF: average O/F ratio (-)

%DEMO-INPUTS%
%{
Demo 1: 
oxMdotHist = linspace(4,1.34,1001);
[totalImp,grainDport,meanOtoF] = getThrustCurve(0:0.01:10,oxMdotHist,1.5,16,1)

Demo 2: 
oxMdotHist = linspace(4,1.34,1501);
[totalImp,grainDport,meanOtoF] = getThrustCurve(0:0.01:15,oxMdotHist,2,11,1)
%}

%INITIALIZATION%
endStep = length(timeHist);
deltaTHist = zeros(1,endStep);
oxFluxHist = zeros(1,endStep);
regRateHist = zeros(1,endStep);
grainDportHist = zeros(1,endStep);
grainDportHist(1) = grainDport;
grainMdotHist = zeros(1,endStep);
OtoFHist = zeros(1,endStep);
chambPressHist = zeros(1,endStep);
thrustHist = zeros(1,endStep);

%MAIN LOOP%
for step = 1:endStep
    if step < endStep
        deltaTHist(step) = timeHist(step + 1) - timeHist(step);
    else
        deltaTHist(step) = timeHist(endStep) - timeHist(endStep - 1);
    end
    [oxFluxHist(step),regRateHist(step),grainMdotHist(step),OtoFHist(step),chambPressHist(step),thrustHist(step)] = getMotorParameters(oxMdotHist(step),grainDportHist(step),grainLength,nozzleArea);
    if step < endStep
    grainDportHist(step + 1) = grainDportHist(step) + 2 * regRateHist(step) * deltaTHist(step);
    end
end

%OUTPUTS%
grainDport = grainDportHist(endStep);
thrustHist(1) = thrustHist(2);
totalImp = sum(thrustHist .* deltaTHist);
meanOtoF = sum(OtoFHist .* deltaTHist) / timeHist(end);

subplot(4,1,1)
plot(timeHist,regRateHist,'r-')
title('Regression Rate (in/sec)')

subplot(4,1,2)
plot(timeHist,grainDportHist,'b-')
title('Grain Port Diameter (in)')

subplot(4,1,3)
plot(timeHist,oxMdotHist,'k-',timeHist,grainMdotHist,'y-')
title('Oxidizer and Fuel Mass Flow Rate (lb/sec)')

subplot(4,1,4)
plot(timeHist,OtoFHist,'g-')
title('O/F (-)')
xlabel('Time (sec)')