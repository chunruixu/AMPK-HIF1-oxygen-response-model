function [] = PlotResult(tout, yout, BeginTime, EndTime, Setting, InsertTime, PlotWindow)

tspan = [0 tout(end)];
[timeL variableNum] = size(yout);

fontsize1 = 18; fontsize2 = 20;

%% plot %%
% O2 = y(1);
% ROS = y(2);
% NADH = y(3);
% NAD = y(4);
% deltaH = y(5);
% ADP = y(6);
% pAMPK =y(7);

% dydt = zeros(14, 1);
% O2 = y(1); 1 对
% pAMPK = y(2); 0.05 对
% ROS = y(3); 0.5 对
% SCAV = y(4);%ROS scavenger 0.1 对
% deltaH = y(5); 0.1 对
% AMP = y(6); 3 
% HIF1a_free = y(7); 1
% HIF1a_AC = y(8); 0.1
% HIF1a_OH = y(9); 0.1
% HIF1 = y(10); 0.1
% NAM = y(11); 1
% NAD = y(12); 3 对 可改为0.5
% NADH = y(13); 0.03 对  可改为0.1
% SIRT1 = y(14); 0.1
%   y0 = [ 1.0000    0.0500    0.5000    0.1000    0.1000    3.0000    1.0000    0.1000    0.1000    0.1000    1.0000    3.0000    0.0300    0.1000]
global Atot NAtot AMPKtot
NAtotyout = zeros(timeL, 1);
NAtotyout = NAtotyout + NAtot;
SIRT1__NAM = NAtotyout - yout(:, 11) - yout(:, 12) - yout(:, 13);
Atotyout = zeros(timeL, 1) + Atot;
ATP = Atotyout - yout(:,6);
%% hyperoxia data
dpHIF1a = [0	 0.097; 12  0.029; 24  0.033];%ref HB10 - H1299
% dpHIF1a_AC = [0	 0.1085; 12  0.0375; 24  0.034];%ref HB10 - H1299 + a549 average
dpHIF1a2 = [0 1; 0.5 0.86; 2 0.78];%RG11
dpSIRT1 = [0 19366; 16 25786];%RG5
dpATP = [0 1; 24 0.99; 48 0.78];%rg29 -- 好像不是人类

%% normal + Hyperoxia %%
ind1 = find(tout>=BeginTime); ind2 = find(tout>=EndTime);
TOUT = tout(ind1:ind2); YOUT = yout(ind1:ind2,:);



figure()%
set(gca,'FontName','Times New Roman','FontSize',16);
colororder({'k','#0000CD'})
   yyaxis left
p2 = line(TOUT, YOUT(:, 6), 'Color','k', 'LineWidth', 2, 'Linestyle', '--')

p3 = line(TOUT, ATP(ind1:ind2), 'Color','k', 'LineWidth', 2, 'Linestyle', '-')

hold on
xa = [InsertTime 45]; %%%%%%%%%%%%%%%%%%%%%%%INSERT ARROW; TUNE HERE
ya = [InsertTime 35];
drawarrow(xa,ya,'textarrow',gca,'Hyperoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
  yyaxis right
  p4 = line(TOUT, YOUT(:,6)./ATP(ind1:ind2), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-');
p9 = line(TOUT, YOUT(:, 14), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-.');

h = legend('ADP', 'ATP','ADP/ATP','free SIRT1','Location', 'NorthEast','fontsize',fontsize1);
xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
ylim([0 1.5]);
box on


figure()%
set(gca,'FontName','Times New Roman','FontSize',16);
colororder({'k','#0000CD'})
   yyaxis left
  pAMPKsim = YOUT(:,2); pAMPKfold = pAMPKsim./pAMPKsim(1);
p4 =line(TOUT, YOUT(:,6)./ATP(ind1:ind2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p9 = line(TOUT, pAMPKfold, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-.');
hold on
xa = [InsertTime 5.5]; %%%%%%%%%%%%%%%%%%%%%%%INSERT ARROW; TUNE HERE
ya = [InsertTime 4.8];
drawarrow(xa,ya,'textarrow',gca,'Hyperoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
  yyaxis right
 p2 = line(TOUT, YOUT(:, 6), 'Color','#0000CD', 'LineWidth', 2, 'Linestyle', '--')
p3 = line(TOUT, ATP(ind1:ind2), 'Color','#0000CD', 'LineWidth', 2, 'Linestyle', ':')
h = legend('ADP/ATP','pAMPK fold','ADP', 'ATP','Location', 'NorthEast','fontsize',fontsize1);
xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
% ylim([0 1.5]);
box on




     figure()%NAD species
 set(gca,'FontName','Times New Roman','FontSize',16);
 colororder({'k','#0000CD'})
   yyaxis left
     p12 = line(TOUT, YOUT(:, 12)./YOUT(:, 13), 'Color',  'k', 'LineWidth', 2, 'Linestyle', '-');

    hold on
xa = [InsertTime 42];
ya = [InsertTime 35];
drawarrow(xa,ya,'textarrow',gca,'Hyperoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
% ylim([0 2.5]);
  yyaxis right
      p5 = line(TOUT, YOUT(:, 13), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '--');
    p5 = line(TOUT, YOUT(:, 11), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', ':');
  p16 = line(TOUT, YOUT(:, 12), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-.');

 h = legend(  'NAD^{+}/NADH','NADH', 'NAM','NAD^{+}', 'Location', 'NorthWest','fontsize',fontsize1);
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
% ylim([0 0.75]);
% ylim([0 2.5]);
box on

  

figure()%SIRT, HIF
set(gca,'FontName','Times New Roman','FontSize',16);
colororder({'k','#0000CD'})
   yyaxis left
 p9 = line(TOUT, YOUT(:, 14), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
%   p15 = line(TOUT, YOUT(:, 9), 'Color', 'k', 'LineWidth', 2, 'Linestyle', ':');

hold on
xa = [InsertTime 0.7]; 
ya = [InsertTime 0.65];
drawarrow(xa,ya,'textarrow',gca,'Hyperoxia','[0.55 0.27 0.07]');
  yyaxis right
 p14 = line(TOUT, YOUT(:, 8), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-');
%    p16 = line(TOUT, YOUT(:, 10), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-.');
 h = legend('SIRT1','HIF1\alpha-AC',  'Location', 'NorthWest','fontsize',fontsize1);
xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
% ylim([0 1.5]);
box on

figure()%SIRT, HIF
set(gca,'FontName','Times New Roman','FontSize',16);
colororder({'k','#0000CD'})
   yyaxis left
   sim = YOUT(:, 14); ini = sim(1);
 p9 = line(TOUT, YOUT(:, 14)./ini, 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
%   p15 = line(TOUT, YOUT(:, 9), 'Color', 'k', 'LineWidth', 2, 'Linestyle', ':');
 sim = YOUT(:, 10); ini = sim(1);
   p16 = line(TOUT, YOUT(:, 10)./ini, 'Color', 'k', 'LineWidth', 2, 'Linestyle', ':');
  

   ylim([0 2.5]);
hold on
xa = [InsertTime 2.5]; 
ya = [InsertTime 2];
drawarrow(xa,ya,'textarrow',gca,'Hyperoxia','[0.55 0.27 0.07]');
  yyaxis right

 p14 = line(TOUT, YOUT(:, 8), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-');

 h = legend('SIRT1 Change Fold', 'HIF1 Change Fold','HIF1\alpha-AC', 'Location', 'NorthWest','fontsize',fontsize1);
xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);

box on


figure()%SIRT, HIF
set(gca,'FontName','Times New Roman','FontSize',16);
colororder({'k','#0000CD'})
   yyaxis left
   p8 = line(TOUT, SIRT1__NAM(ind1:ind2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', ':');
hold on
% xa = [InsertTime 0.7]; %%%%%%%%%%%%%%%%%%%%%%%INSERT ARROW; TUNE HERE
% ya = [InsertTime 0.65];
% drawarrow(xa,ya,'textarrow',gca,'Hyperoxia','[0.55 0.27 0.07]');
  yyaxis right
 p14 = line(TOUT, YOUT(:, 10), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-');
%    p16 = line(TOUT, YOUT(:, 10), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-.');
 h = legend('SIRT1:NAM','HIF1',  'Location', 'NorthWest','fontsize',fontsize1);
xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
% ylim([0 1.5]);
box on
% figure()
%  set(gca,'FontName','Times New Roman','FontSize',fontsize2);
% p9 = line(TOUT, YOUT(:, 14), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
% p8 = line(TOUT, SIRT1__NAM(ind1:ind2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', ':');
% 
% %  p14 = line(TOUT, YOUT(:, 8), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-');
% %   p15 = line(TOUT, YOUT(:, 9), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '--');
%  p14 = line(TOUT, YOUT(:, 8), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-.');
%   p15 = line(TOUT, YOUT(:, 9), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
%    hold on
% xa = [InsertTime 0.9];
% ya = [InsertTime 0.75];
% drawarrow(xa,ya,'textarrow',gca,'Hypoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
% h = legend('free SIRT1', 'SIRT1:NAM','HIF1\alpha-AC', 'HIF1\alpha-OH', 'Location', 'NorthWest','fontsize',fontsize1);
%      xlabel('Time (h)','fontsize',fontsize1);
% box on

     
%% specific setting plot %%
ind11 = find(tout>=InsertTime); ind22 = find(tout>=PlotWindow);
YOUT2 = yout(ind11:ind22,:); TOUT2 = tout(ind11:ind22); TOUT2 = TOUT2-TOUT2(1);
if strcmp(Setting,'hyperoxia1')

figure()%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'FontName','Times New Roman','FontSize',16);
HIF1ACsim = YOUT2(:,8);
[score,scalary] = findSquares(TOUT2,HIF1ACsim,dpHIF1a);
[score2,scalary2] = findSquares(TOUT2,HIF1ACsim,dpHIF1a2);
line(TOUT2, HIF1ACsim, 'Color', 'k', 'LineWidth', 2);
hold on;
scatter(dpHIF1a(:,1),dpHIF1a(:,2).*scalary, 100,'ro','LineWidth', 2);
hold on;
scatter(dpHIF1a2(:,1),dpHIF1a2(:,2).*scalary2, 100,'r^','LineWidth', 2);
xlabel('Time after hyperxia (h)','fontsize',fontsize1); ylabel('HIF1\alpha-AC','fontsize',fontsize1);
 h = legend( 'simulation', 'Liu et al.','Haugen et al.','Location', 'NorthWest','fontsize',fontsize1);
 xlim([0 50]);
 box on;

 figure()%%%%%%%%%%%%%%%%%%%%%%%%%%%%SIRT1
 set(gca,'FontName','Times New Roman','FontSize',16);
 SIRT1sim = YOUT2(:,14);
[score,scalary] = findSquares(TOUT2,SIRT1sim,dpSIRT1);
     line(TOUT2, SIRT1sim, 'Color', 'k', 'LineWidth', 2);
     hold on;
     scatter(dpSIRT1(:,1),dpSIRT1(:,2).*scalary, 100,'ro','LineWidth', 2);
      xlabel('Time after hyperxia (h)','fontsize',fontsize1); ylabel('SIRT1','fontsize',fontsize1);
   h = legend( 'simulation', 'Podgorski et al.','Location', 'NorthWest','fontsize',fontsize1);
xlim([0 50]);
box on;

 figure()%%%%%%%%%%%%%%%%%%%%%%%%%%%%ATP
     set(gca,'FontName','Times New Roman','FontSize',fontsize2);
     ATPsim = Atot-YOUT2(:,6);
     line(TOUT2, ATPsim, 'Color', 'k', 'LineWidth', 2);
     hold on;
     [score,scalary] = findSquares(TOUT2,ATPsim,dpATP);
     scatter(dpATP(:,1),dpATP(:,2).*scalary,100,'r','o','LineWidth', 2);
%      hold on;
%         [score,scalary] = findSquares(TOUT2,ATPsim,dpATP4);
%      scatter(dpATP4(:,1),dpATP4(:,2).*scalary, 100, 'r','o', 'LineWidth', 2);
     xlabel('Time after hypoxia (h)','fontsize',fontsize1);ylabel('ATP','fontsize',fontsize1);
h = legend( 'simulation', 'Gille et al.', 'Location', 'North','fontsize',fontsize1);%ref76;86
% h = legend( 'simulation', 'Casey et al. 2002', 'Location', 'North','fontsize',fontsize1);%ref76;86
box on



     
        figure()
set(gca,'FontName','Times New Roman','FontSize',16);
   NADH = YOUT2(:,13); 
    line(TOUT2, NADH, 'Color', 'r', 'LineWidth', 2);
%      h = legend( 'sim', 'data', 'Location', 'North');
     xlabel('Time (h)'); ylabel('NADH');
     
     
     figure()%pAMPK fold
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
pAMPKsim = YOUT2(:,2); pAMPKfold = pAMPKsim./pAMPKsim(1);
line(TOUT2, pAMPKfold, 'Color', 'k', 'LineWidth', 2);
hold on;
% scatter(dppAMPK(:,1),dppAMPK(:,2),50,'ko','filled');
% scatter(dppAMPK1(:,1),dppAMPK1(:,2),100,'r','o','filled');
% scatter(dppAMPK2(:,1),dppAMPK2(:,2),50,'r','o','filled');
xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('pAMPK fold','fontsize',fontsize1);
%  h = legend( 'simulation', 'Wu et al. and Chen et al.', 'Location', 'North','fontsize',fontsize1);%r106-34   
%  ylim([1 3.2]);
box on 

end




