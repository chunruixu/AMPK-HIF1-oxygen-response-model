clear all
%multiple rounds of parameter search
max_run_time = 5*3600; %set the maximum time (s)
%input an intial set of parameters
[param0] = getParam(4); 
y0 = [ 1.0000    0.0500    0.5000    0.1000    0.1000    3.0000    1.0000    0.1000    0.1000    0.1000    1.0000    3.0000    0.0500    0.1000]; 

%% LOAD params and initial y0
load('Param_Collection/2025-06-05 22.25.12 std_0 Svalue_5.8787.mat')%%%一种同上 AMPK和SIRT差一点
% load('Param_Collection/2025-06-08 21.31.07 std_0 Svalue_6.6819.mat')
load('Param_Collection/2025-06-09 10.16.55 std_0 Svalue_7.4579.mat')
load('Param_Collection/2025-06-09 12.49.17 std_0.01 Svalue_3.0895.mat')
load('y0update3.mat')
% load('Param_Collection/2025-06-09 16.54.11 std_0 Svalue_7.2009.mat')



param0 = Pbest;%!!!!


y0List = [1:14]; 
%determine conditions you want to consider in parameter estimation
% settingList = ["normoxia","HighpAMPK1","hypoxia1"];
settingList = ["normoxia","hypoxia1"];
% settingList = ["normoxia","HighpAMPK1"];
% settingList = ["normoxia"];
% settingList = ["normoxia","hypoxia1","hyperoxia1"];

timeBegin = 700; timeEnd = 1000;%to check normal condition
[cost, costStorage] = getCost(y0, param0, timeBegin, timeEnd, settingList);%compute the cost for the initial parameters

% [Pbest, y0New, Sbest]=ParameterEstimationGA(y0,y0List, stdv_y, param0, settingList, iter1, PopS, timeBegin,TimeOfRun)  
% tic


% [y0] = randomPar(y0,y0List,0.01);
std = 0.02;
for i=1:100
    tic %record time
%genetic algorithm
    [Pbest, y0New, Sbest]=ParameterEstimationGA(y0, y0List, std, param0, settingList, 200, 50, timeBegin,timeEnd) ;%maxG, pSize
    %keep y0New as the new intial values for the next round of
    %parameterizatioin if the cost is less/equal to a threshold
    if Sbest == cost
        std = std + 0.02;
    else
        std = max((std-0.01),0);
    end
    if Sbest<=cost*1.2 
        y0 = y0New; param0 = Pbest; cost = Sbest;
    end
    toc%check time
    if toc > max_run_time 
        error("Time Out !!!"); 
    end 
    
end
