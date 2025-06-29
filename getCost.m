% function cost = getCost(yout, tout, timeBegin, timeEnd, setting)
function [COST costStorage] = getCost(y0, param, timeBegin, timeEnd, settingList)
% function cost = getCost(y0, param, timeBegin, timeEnd, settingList)
% at normal condition - first step - ensure the system is relatively reasonable 
%O2 1-10%
 global Atot NAtot
param0 = param;
sim = length(settingList);%number of conditions considered
costStorage = zeros(sim,1);
for i=1:sim%%parfor
    setting = settingList(i);
    if strcmp(setting,'normoxia') %
        [tout, yout] = Sim(timeEnd,param0, y0);
        [AverageO2, ss3, ~] = getAverage(yout(:, 1), tout, timeBegin, timeEnd);
        cost1 = (max(abs(AverageO2 - 7.75) - 7.25, 0))^2;
        [AverageNADRatio, ~, ~] = getAverage(yout(:,12)./yout(:,13), tout,timeBegin, timeEnd);
        cost2 = (max(abs(AverageNADRatio - 350.5) - 349.5, 0))^2; %1~700
        [AverageSumNAD, ~, ~] = getAverage(yout(:,12)+yout(:,13), tout,timeBegin, timeEnd);
        cost3 = (min((AverageSumNAD/NAtot - 0.2), 0))^2;
        Atotyout = zeros(length(tout), 1) + Atot;
        ATP = Atotyout - yout(:,6);
        [AverageATPRatio, ss1, ~] = getAverage(ATP./yout(:,6), tout,timeBegin, timeEnd);
         cost4 = (max(abs(AverageATPRatio - 10.5) - 9.5, 0))^2;
        [AverageHIF1a_AC, ss2, ~] = getAverage(yout(:,8), tout,timeBegin, timeEnd);
        
        if ss3 == 0
            cost5 = 1;
        else
            cost5 = 0;
        end
        if ss2 == 0
            cost6 =1;
        else
            cost6 = 0;
        end
         cost = (cost1 + cost2 + cost4+cost3)*10 +100*cost5;  %normal---RECORD
    elseif strcmp(setting,'hypoxia1')
        [tout, yout] = Sim(timeEnd,param0, y0);
         y0 = yout(end,:);%run normal sim and get the y0 at steady-state
         [param] = getMutantParam(param0,'hypoxia1');
         [tout2, yout2] = Sim(timeEnd,param, y0);
        O2level1 = getAverage(yout(:,1), tout,timeBegin, tout(end));
        O2level2 = getAverage(yout2(:,1),tout2,0,tout2(end));
     cost7 = (max((abs((O2level2/O2level1)-0.075)-0.025),0))^2;%O2 decrease to 1/20~1/10; 1%~2%
     
     pAMPKsim = yout2(:,2); 
%         pAMPKfold = pAMPKsim(end)/pAMPKsim(1);
        pAMPKfold = pAMPKsim./pAMPKsim(1);%0905 - ref66
cost8 = (min((pAMPKfold(end)-2.8),0))^2;%more than 2fold %0428 
 cost8 = (max(abs(pAMPKfold(end) - 3.75) - 1.25, 0))^2; %限制一下最大值；没有参考文献
dppAMPK = [0 1; 2 1.5; 4 2;0 1; 2 1.5; 4 2; 12 3]; %REF66 -tuned by ref68%0428弱化最后一个点
        cost9 = findSquares2(tout2,pAMPKfold,dppAMPK);
%         cost13 = 300*findSquares(tout,pAMPKfold,dppAMPK);
%    [~, ss4, ~] = getAverage(yout2(:, 2), tout2, 0,tout2(end));
%         if ss4 == 0
%             cost15 =1;
%         else
%             cost15 = 0;
%         end

   dpNAD = [0 1; 8 0.8; 16 0.5];%R106REF34  -fig5
   dpNAD = [0 8.47; 4 7.97; 12 6.96; 24 5.1]; %R106ref34 fig7
  NAD =  yout2(:,12);
  cost10 = findSquares(tout2,NAD,dpNAD);
   NADsim = yout2(:,12); 
          NADfold = NADsim./NADsim(1);%
          if NADfold(end)<=1
              cost14 = 0;
          else
              cost14 = 1;
          end
%HIF%    HIF1a = yout2(:,1)+yout2(:,2)
    HIF1a_AC = yout2(:,8);
% HIF1a_tot = yout2(:,7)+ yout2(:,8)+ yout2(:,9);
    dpHIF1a_AC = [0	 0.06; 4  1.85; 8  2.3; 16  0.5];% 16  0.5; 16  0.5];%ref108
%     dpHIF1a_AC = [0	 1; 4  7; 4  7; 8  10];%24,9.1; 48 2.9];%ref108:fold; HEK293 mRNA  %parameterization
    cost11 = findSquares(tout2,HIF1a_AC,dpHIF1a_AC);
% cost11 = findSquares(tout2,HIF1a_AC,dpHIF1a3);
    SIRT1 = yout2(:,14); %%%%%%%%%%%%%%%%加SIRT:NAM
      dpSIRT1 = [0 0.520; 2 0.754; 4 1.254; 8 0.932];
      dpSIRT1 = [0 1.084086763;24 0.173202546; 8 0.7];%R3REF34
      cost12 = findSquares(tout2,SIRT1,dpSIRT1);   
%%13 atp
Atotyout = zeros(length(tout2), 1) + Atot;
ATP = Atotyout - yout2(:,6);
dpATP =[0 16.63; 0.5 14.54; 0.75 12.68; 1 10.81];%r76 
dpATP2 = [0 24.77; 0.5 17; 1 11.61; 1.5 1.1];%R78  
dpATP3 = [0 17.99;0.25	13.49;0.5	4.63;0.75	3.41;1	2.41];%R82
dpATP4 = [0	25.76;1.5	18.28;3	14.27;5	8.05];%R86
cost13 = findSquares(tout2,ATP,dpATP);



  cost = 100*cost7 +20*cost8+50*cost9 +  20*cost11 +20*cost12+ 200*cost13 +5*cost14;%9:pAMPK; 11:HIF;12:SIRT  + 5*cost10(10 NAD+)

    elseif strcmp(setting,'hyperoxia1')
     %%SIRT
     dpSIRT1 = [0 19366; 16 25786];%RG5

     SIRT1 = yout2(:,14); 
      cost15 = findSquares(tout2,SIRT1,dpSIRT1);   
%      SIRT1fold = SIRT1./SIRT1(1);
%      if SIRT1fold > 1;
%      cost15 = 0;
%      else
%          cost15 =1;%RG3-5;与RG1不同
%      end
      %%HIF
     HIF1a_AC = yout2(:,8);
     HIF1a_ACfold = HIF1a_AC./HIF1a_AC(1);
     if HIF1a_ACfold < 1;
     cost16 = 0;
     else
         cost16 =1;%RG1,rg18,rg20
     end
     cost = 100*cost15+5*cost16;
    end
    costStorage(i) = cost;
end
COST = sum(costStorage);
% cost = cost1 + cost2*1000 + cost4*1000;
% cost = cost4;

