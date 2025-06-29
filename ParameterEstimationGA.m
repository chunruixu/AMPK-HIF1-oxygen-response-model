function [Pbest, y0New, Sbest]=ParameterEstimationGA(y0,y0List, stdv_y, param0, settingList, iter1, PopS, timeBegin,TimeOfRun)    
%% [Pbest, ~]=ParameterEstimationGA(y0,param, settingList, 100, 10, 50,150)   
%  settingList = ["normoxia","hypoxia"];
 %paramList = [1:24]; which param is tuned; stdv: see 'randomPar'
 global AMPKtot NAtot
 %define LB and UB
% param = [1 1 0.1 1 1 10 0.5 30 10 1 0.5 0.1 1 0.1 0.5 0.5 0.5 1 0.5 10 0.2 0.03 5 0.01];
%param = [ 1.0000    0.0500    0.5000    0.1000    0.1000    3.0000    1.0000    0.1000    0.1000    0.1000    1.0000    3.0000    0.0500    0.1000]; 
LB = 0.05.*param0;
UB = 20.*param0;
LB(10) = 0.3; UB(10) = 30; % %R43 ~3
LB(50) = 0; UB(50) = 0.9;% alpha
LB(4) = 0.01; UB(4) = 1; %Jm_etco2 %R43 ~table3  0.1
% UB(8) = param0(8);%ROS-NADH通过NOX
%%%%%%%%
% LB = 0.1.*param0;
% UB = 10.*param0;
% UB(26) = param0(26);
% LB(25) = param0(25);


% LB(9) = 0.1*AMPKtot; UB(9) = 10*AMPKtot;%Ji_ampknox
% LB(20) = 0.005; UB(20) = 100;%Ja_adpampk
% LB(6) = 0.01; UB(6) = 10000;%Ji_nado2
% % NAtot = param0(34);
% LB(21) = 0.1*NAtot; UB(21) = 10*NAtot; %Ja_nadhca


% LB(34) = 1; UB(34)=10;
% LB(34) = 4; UB(34)=4;
%k_o2input, k_etc_f, k_etc_b, k_nado2, Ji_nado2, k_nox_f, ...
% Ji_ampknox, k_nox_b, kd_ros, k_nam_f, k_nam_b, 
% k_a_f, k_a_b, k_ampk_phos1, k_ampk_phosLKB, k_ampk_phosCaM,
% Ja_nadhca, k_ampk_unphos, Jm_etco2, Jm_noxo2, alpha,
% Ja_adpampk, k_nad_f, k_nad_b

function cost = FUNGA(param)
    try
    [cost, costStorage] = getCost(y0, param, timeBegin, TimeOfRun, settingList);
    catch
        cost = inf;
    end
end

options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf,'InitialPopulation', param0,'MaxGenerations',iter1,'PopulationSize', PopS);


[newparam, fval] = ga(@FUNGA,length(param0),[],[],[],[],LB,UB,[],options)

Pbest = newparam; 
Sbest = fval;
%generate a new y0
if stdv_y ~= 0
[y0New] = randomPar(y0,y0List,stdv_y);
else
    y0New = y0;
end


datetim=datetime('now');
DateString = char(datetim);
DateString(DateString==':')='.';

%
foldername = 'Param_Collection';
% foldername = 'Param_Collection1';
filename=[foldername,'/',DateString, ' ','std_',num2str(stdv_y),' Svalue_', num2str(round(Sbest,4)),'.mat'];
save(filename,'Pbest','Sbest', 'y0', 'settingList')
% 
% save()
end