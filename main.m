clear all %clear temporary variables
close all ;%close all open figures
global Atot NAtot AMPKtot 
%set some global variables for [ATP+AMP], [NAD+NADH+NAM], and [AMPK+pAMPK]
Atot = 30;% 30e3; %[R31] TABLE3; 30mM
NAtot = 2*2;%15; %[R31] TABLE3 2mM;   %[R13] physiological intracellular NADH: 1-10uM=0.001-0.01mM   %NAM in mammalian tissues - R7 - 11~400uM
AMPKtot = 0.1;%3.2e-4mM;%20ppm - pax database, integrated whole organism
% NAtot = param(34);

%set a preliminary initial values for each variable
y0 = [ 1.0000    0.0500    0.5000    0.1000    0.1000    3.0000    1.0000    0.1000    0.1000    0.1000    1.0000    3.0000    0.0500    0.1000]; 
%initial conditions for O2 pAMPK ROS SCAV deltaH AMP HIF1a_free HIF1a_AC HIF1a_OH HIF1 NAM NAD NADH SIRT1 
[param] = getParam(); %load initial parameters  

%% OR load specific parameters and initial conditions obtained from parameterizatioin as follows:
load('Param_Collection/2025-06-13 11.28.10 std_0 Svalue_4.9273.mat')
load('y0_4.9273.mat')
param = Pbest; 
y0=y0;


%% sim at different stages
% TimeOFChange = 100; %run time
TimeOfRun_norm = 24;%150;%>=50
TimeOfRun_hyp = 72;
y0_1 = y0; %the initial condition of normal situation
[tout1, yout1] = Sim(TimeOfRun_norm,param, y0_1); %simulation of normal situation

y0_2 = yout1(end,:); %the initial condition of hypoxic situation or others
%% switch different settings - hypoxia or hyperoxia
mutant = 'hypoxia1';%'HighpAMPK1';
% mutant = 'hyperoxia1';%
% mutant = 'HighpAMPK1';
% mutant = 'HighSIRT1';
% mutant = 'NMN';
% mutant ='ADPATP';
% mutant = 'AMPKinhibitor&hypoxia';
% mutant = 'HighpAMPKALL';

[param] = getMutantParam(param,mutant);
[tout2, yout2] = Sim(TimeOfRun_hyp,param, y0_2);
 yout = [yout1(1:(end-1),:); yout2];
 tout = [tout1(1:(end-1)); tout2+TimeOfRun_norm];

 %CHECK ALL DYNAMICS
%   PlotResult(tout, yout, 0, TimeOfRun_norm+TimeOfRun_hyp)

if strcmp(mutant,'hypoxia1')
 PlotResult2(tout, yout, 1, TimeOfRun_norm+TimeOfRun_hyp, mutant, TimeOfRun_norm, TimeOfRun_norm+48)%timecourese %lvl, lvltime1, lvltime2  'hypoxia1'
end

if strcmp(mutant,'hyperoxia1')
 PlotResult3(tout, yout, 1, TimeOfRun_norm+TimeOfRun_hyp, mutant, TimeOfRun_norm, TimeOfRun_norm+48)%timecourese %lvl, lvltime1, lvltime2  'hypoxia1'
end

if strcmp(mutant,'HighpAMPK1')||strcmp(mutant,'HighSIRT1')||strcmp(mutant,'HighpAMPKALL')
    HIFACfold = yout(end,8)/yout1(end,8)
    pAMPKfold = yout(end,2)/yout1(end,2)
    pAMPK = yout(end,2)
    HIFAC = yout(end,8)
    HIF1 = yout(end,10)
    NADratio = yout(end,12)/yout(end,13)
    SIRT1 = yout(end,14)
    SIRT1fold = yout(end,14)/yout1(end,14)
 PlotResult4(tout, yout, 10, TimeOfRun_norm+TimeOfRun_hyp, mutant, TimeOfRun_norm, TimeOfRun_norm+50)%timecourese %lvl, lvltime1, lvltime2  'hypoxia1'
end

%%
