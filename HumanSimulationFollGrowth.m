%set specifications for simulations a call calculation and solving of the ODE
%needed:'Model28_Parameter.dat', 'InitialValues.txt'

function HumanSimulationFollGrowth
%
%-----------------------------------------------------------------------
%
runnum        = 1;
%save simulation results
ShowStuff     = 1;
SaveStuff     = 0;
DirStuff      = 'Results_LutStim/'; 
%select type of simulation 
NormalCycle   = 0;
LutStim       = 0;
LateFollPhase = 1; 
%
%-----------------------------------------------------------------------
% 
for i = 1:runnum
%
%-----------------------------------------------------------------------
%   
%integration time beginning and end
    tb = 0;
    te = 200;
%
%-----------------------------------------------------------------------
% 
% several parameter vectors (function see ReadMe)
para=[];
    para(1) = 0;                %ODE function called to test(0) or not (1)
    para(2) = 54;               %number of extra-equations 
    para(3) = 0;
    para(4) = 0;           
    para(5) = 0;   
    para(6) = 0;    
    para(7) = 6/10;             %mean for FSH Sensitivity 
    para(8) = 0.55*para(7);     %std.deviation for FSH Sensitivity %0.55
    para(9) = 25;               %threshold LH concentration for ovulation    
    para(10) = 5;               %big but not ovulated follicle livetime
    para=para';

%parameter for growth equation
    v       = 2;                      %fractal dimension 
    gamma   = 0.035/2;                %growth rate %0.035/2
    xi      = 25;                     %max. diameter of follicles 
    mu      = 1;                      %proportion of self harm  
    k       = 0.065/(xi^v);           %strength of competition %0.065
    rho     = 0.01;                   %rate of decline
    Folmax  = 22;                     %min. ovulation size %20
paraOde = [ v; gamma; xi; mu; k; rho; Folmax];

%parameters for poisson distribution
lambda          = 20/14;           %#Follikels/days
intervallPoi    = 0.25;            %intervall per day in which follicles appear
paraPoi = [lambda; intervallPoi];
  
%load parameter for ODE
file = 'Parameter.dat';
Par=importdata(file,' ');
%
%-----------------------------------------------------------------------
%
%load initial values
file = 'InitialValues.txt';
delimiterIn=';';
headerlinesIn=0;
yInitial=importdata(file,delimiterIn,headerlinesIn);
%
%-----------------------------------------------------------------------
%
%initial follicles
y0Foll = 4;                                                 
StartValues = [y0Foll; yInitial]';

[FSHVec, StartVec] = CreateFollicles(para,paraPoi,tb,te);
%
%-----------------------------------------------------------------------
%
%Normal Cycle
%
%-----------------------------------------------------------------------
%
if NormalCycle
    Simulation_NormalCycle(para,paraPoi,paraOde,Par,tb,te,StartValues,StartVec,FSHVec,ShowStuff,SaveStuff,DirStuff);
end
%
%-----------------------------------------------------------------------
%
%Luteal Phase Stimulation: FSH/LH administartion (Menopur)
%
%-----------------------------------------------------------------------
%
if (LutStim)
    Par(599)=200;          %start of dosing - fiktive Zeitpunkte, werde in der Simulation gesetzt
    Par(598)=Par(599)+10;   %end time of dosing       
    Par(579)=13.387/2.65;   %D FSH
    Par(578)=9.87;          %beta FSH
    Par(577)=0.42;          %clearance rate FSH
    Par(576)=2.14;          %D LH
    Par(575)=6.04;          %beta LH 
    Par(574)=3.199;         %clearance rate LH
    Par(596) = 0;
    Simulate_LutStim(para,paraPoi,paraOde,Par,tb,te,StartValues,StartVec,FSHVec,ShowStuff,SaveStuff,DirStuff);
end
%
%-----------------------------------------------------------------------
%
%Late Folliculare Phase Phase Stimulation: FSH/LH administartion (Menopur)
%
%-----------------------------------------------------------------------
%
if (LateFollPhase)
    Par(599)=200;          %start of dosing - fiktive Zeitpunkte, werde in der Simulation gesetzt
    Par(598)=Par(599)+10;   %end time of dosing       
    Par(579)=13.387/4;      %D FSH
    Par(578)=9.87;          %beta FSH
    Par(577)=0.42;          %clearance rate FSH
    Par(576)=2.14;          %D LH
    Par(575)=6.04;          %beta LH 
    Par(574)=3.199;         %clearance rate LH
    Par(596) = 0;
    Simulate_LateFollPhase(para,paraPoi,paraOde,Par,tb,te,StartValues,StartVec,FSHVec,ShowStuff,SaveStuff,DirStuff);
end
%
%-----------------------------------------------------------------------
%
end