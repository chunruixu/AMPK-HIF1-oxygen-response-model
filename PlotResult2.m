function [] = PlotResult(tout, yout, BeginTime, EndTime, Setting, InsertTime, PlotWindow)

tspan = [0 tout(end)];
[timeL variableNum] = size(yout);

fontsize1 = 22; fontsize2 = 18;

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
%% data
 %r6-ref44-fig4 data, C2C12 cells
NAMdata = [42.5 11.3]; NADratiodata = [6 17];
 %r-m2-ref170-fig4  rat heart
 AMPratiodata1 = [0.04 0.11 0.15]; pAMPKdata1 = [107.53 672.04 685.48];
 %r-m2-ref236 H4IIE cells
 ADPratiodata2 = [0.035 0.052 0.056 0.068 0.090 0.093 0.104]; pAMPKdata2 = [0.139 0.149 0.161 0.241 0.507 0.451 0.396];
 %integrated data from the two above paper
%   ADPratiodata3 = [0.070 0.104 0.112 0.136 0.181]; pAMPKdata3 = [1 1.271 1.157 1.734 3.651];
ADPratiodata3 = [0.070 0.104 0.112 0.136 0.181]; pAMPKdata3 = [0.139 0.149 0.161 0.241 0.507];
  %%temporal data
dppAMPK = [0 1; 2 1.5; 4 2; 12 3]; %REF66 -tuned by ref68, hypoxia
dppAMPK1 = [0 1; 2 1.5; 4 2]; %larger weight
dppAMPK2 = [12 3]; %less weight
dpNADratio = [0 2.88; 2 3.68; 4 5.29; 8 5.18]; %ref6-ref45, part, AMPK activator
%    dpNAD1 = [0 0.315; 0.17 0.331; 1 0.442; 2 0.496];%R70
   dpNAD = [0 1; 8 0.8; 16 0.5];%R106REF34
   dpNADH = [0	1; 0.5	0.5; 1	1; 1.2	1.2; 1.5	0.9; 1.7	0.7; 2	0.6; 2.2	0.5; 2.5	0.4; 2.7	0.26;3	0.15];%NADH1
dpNADH2 = [0.083333333 0.52; 0.166666667 0.97; 0.25 1.81; 0.333333333 1.84; 0.416666667 1.82];%R69-参考意义稍微有点小，因为常氧咋也变化啊
dpATP =[0 16.63; 0.5 14.54; 0.75 12.68; 1 10.81];%r76 
dpATP2 = [0 24.77; 0.5 17; 1 11.61; 1.5 1.1];%R78  
dpATP3 = [0 17.99;0.25	13.49;0.5	4.63;0.75	3.41;1	2.41];%R82
dpATP4 = [0	25.76;1.5	18.28;3	14.27;5	8.05];%R86
% dpHIF1a = [0	 0.06; 4  1.85; 8  2.3;12,2.25];%ref108: 綜合幾個fig
dpHIF1a_1 = [0	 1; 4  7; 8  10;24,9.1;48 2.9];%ref108:fold; HEK293 mRNA 

% dpHIF1a3 = [0	 0.06; 4  1.85; 8  2.3; 16  0.5];% 16  0.5; 16  0.5];%ref108: 綜合幾個fig

dpHIF1a2 = [0 1; 0.5 1.74; 2 3.1];%rg11 fold
% dpHIF1a3 = [0 1; 4 12.27400083; 12 14.43031387; 16 11.28902144];%ref108 HIF1a MCF10A

dpHIF1OHP = [0 1; 8 0.841397395; 24 0.973654254; 48 0.961833283];%r133






dpSIRT1_1 = [0 0.520; 2 0.754; 4 1.254];%;8,0.93;16,0.53];
    dpSIRT1_1 = [0 1.084086763;8 0.698867712; 12 0.769845556; 24 0.173202546];%R3REF34
dpSIRT1_2 = [8,0.93;16,0.53];
dpSIRT3 = [0 1; 6 0.65; 9 0.51; 12 0.43];%R134

%% normal + hypoxia %%
ind1 = find(tout>=BeginTime); ind2 = find(tout>=EndTime);
TOUT = tout(ind1:ind2); YOUT = yout(ind1:ind2,:);

figure()%ATP
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
%  colororder({'k','[0.09412 0.4549 0.80392]'})
% yyaxis left
p2 = line(TOUT, YOUT(:, 6), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--')

p4 = line(TOUT, YOUT(:,6)./ATP(ind1:ind2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');

p3 = line(TOUT, ATP(ind1:ind2), 'Color','k', 'LineWidth', 2, 'Linestyle', ':')
hold on
xa = [InsertTime 35]; %%%%%%%%%%%%%%%%%%%%%%%INSERT ARROW; TUNE HERE
ya = [InsertTime 28];
drawarrow(xa,ya,'textarrow',gca,'Hypoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
% yyaxis right
% p3 = line(TOUT, ATP(ind1:ind2), 'Color','[0.09412 0.4549 0.80392]', 'LineWidth', 2, 'Linestyle', '-')
h = legend('ADP','ADP/ATP','ATP', 'Location', 'NorthWest','fontsize',fontsize1);
xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
% ylim([0.1 0.2]);
box on


figure()%NAD species
 set(gca,'FontName','Times New Roman','FontSize',fontsize2);
 colororder({'k','#0000CD'})
   yyaxis left
    p12 = line(TOUT, YOUT(:, 12)./YOUT(:, 13), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
     p5 = line(TOUT, YOUT(:, 13), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
    hold on
xa = [InsertTime 10];
ya = [InsertTime 8.3];
drawarrow(xa,ya,'textarrow',gca,'Hypoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
% ylim([0.13 0.27]);
  yyaxis right
  p16 = line(TOUT, YOUT(:, 12), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '--');
      
%     p5 = line(TOUT, YOUT(:, 11), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', ':');
%  h = legend( 'NAD^{+}/NADH', 'NADH','NAD^{+}', 'NAM', 'Location', 'NorthWest','fontsize',fontsize1);
 h = legend( 'NAD^{+}/NADH', 'NADH','NAD^{+}', 'Location', 'NorthWest','fontsize',fontsize1);
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
% ylim([0 0.75]);
box on


figure()%SIRT species
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
p8 = line(TOUT, SIRT1__NAM(ind1:ind2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', ':');
p9 = line(TOUT, YOUT(:, 14), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
p10 = line(TOUT, YOUT(:, 14)+SIRT1__NAM(ind1:ind2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
hold on
xa = [InsertTime 0.65];
ya = [InsertTime 0.55];
drawarrow(xa,ya,'textarrow',gca,'Hypoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
h = legend( 'SIRT1:NAM', 'free SIRT1', 'total SIRT1','Location', 'NorthWest','fontsize',fontsize1);
xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
box on
% ylim([0 4]);
% 




 figure()%HIF1a_free ;% HIF1a_AC ;% HIF1a_OH;% HIF1 ;
 set(gca,'FontName','Times New Roman','FontSize',fontsize2);
 colororder({'k','#0000CD'})
   yyaxis left
%    p13 = line(TOUT, YOUT(:, 7), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');
%  p13 = line(TOUT, YOUT(:, 7)+YOUT(:, 8)+YOUT(:, 9), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '-');

        p14 = line(TOUT, YOUT(:, 8), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
 p15 = line(TOUT, YOUT(:, 9), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
    hold on
xa = [InsertTime 0.7];
ya = [InsertTime 0.6];
drawarrow(xa,ya,'textarrow',gca,'Hypoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
% ylim([0.24 0.27]);
  yyaxis right
  p16 = line(TOUT, YOUT(:, 10), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '--');
%   h = legend( 'HIF1\alpha-free', 'HIF1\alpha-OH', 'HIF1\alpha-AC', 'HIF1', 'Location', 'NorthWest','fontsize',fontsize1);
 h = legend('HIF1\alpha-AC', 'HIF1\alpha-OH', 'HIF1', 'Location', 'NorthWest','fontsize',fontsize1);
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
% ylim([0 0.75]);
% ylim([0 0.05]);
box on

figure()
 set(gca,'FontName','Times New Roman','FontSize',fontsize2);
p9 = line(TOUT, YOUT(:, 14), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-');
p8 = line(TOUT, SIRT1__NAM(ind1:ind2), 'Color', 'k', 'LineWidth', 2, 'Linestyle', ':');

%  p14 = line(TOUT, YOUT(:, 8), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '-');
%   p15 = line(TOUT, YOUT(:, 9), 'Color', '#0000CD', 'LineWidth', 2, 'Linestyle', '--');
 p14 = line(TOUT, YOUT(:, 8), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '-.');
  p15 = line(TOUT, YOUT(:, 9), 'Color', 'k', 'LineWidth', 2, 'Linestyle', '--');
   hold on
xa = [InsertTime 0.9];
ya = [InsertTime 0.75];
drawarrow(xa,ya,'textarrow',gca,'Hypoxia','[0.55 0.27 0.07]');%drawarrow(startpoint,endpoint,'String','color',linewidth);
h = legend('free SIRT1', 'SIRT1:NAM','HIF1\alpha-AC', 'HIF1\alpha-OH', 'Location', 'NorthWest','fontsize',fontsize1);
     xlabel('Time (h)','fontsize',fontsize1);
box on


%% specific setting plot %%
ind11 = find(tout>=InsertTime); ind22 = find(tout>=PlotWindow);
YOUT2 = yout(ind11:ind22,:); TOUT2 = tout(ind11:ind22); TOUT2 = TOUT2-TOUT2(1);
if strcmp(Setting,'hypoxia1')

figure()%pAMPK fold
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
pAMPKsim = YOUT2(:,2); pAMPKfold = pAMPKsim./pAMPKsim(1);
line(TOUT2, pAMPKfold, 'Color', 'k', 'LineWidth', 2);
hold on;
% scatter(dppAMPK(:,1),dppAMPK(:,2),50,'ko','filled');
scatter(dppAMPK1(:,1),dppAMPK1(:,2),100,'r','o','filled');
scatter(dppAMPK2(:,1),dppAMPK2(:,2),50,'r','o','filled');
xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('pAMPK fold','fontsize',fontsize1);
 h = legend( 'simulation', 'Wu et al. and Chen et al.', 'Location', 'North','fontsize',fontsize1);%r106-34   
 ylim([1 3.2]);
box on 

figure()
set(gca,'FontName','Times New Roman','FontSize',16);
     NADsim = YOUT2(:,12);
     [score,scalary] = findSquares(TOUT2,NADsim,dpNAD);
     line(TOUT2, NADsim, 'Color', 'k', 'LineWidth', 2);
    hold on;
     scatter(dpNAD(:,1),dpNAD(:,2).*scalary, 50,'r','o', 'LineWidth', 2);
      xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('NAD^{+}','fontsize',fontsize1);
      h = legend( 'simulation', 'Lim et al. 2010', 'Location', 'North','fontsize',fontsize1);%r106-34
box on
%      
     figure()%%%%%%%%%%%%%%%%%%%%%%%%%%%NADH
     set(gca,'FontName','Times New Roman','FontSize',16);
% %      ind31 = find(tout>=InsertTime); ind32 = find(tout>=(InsertTime+0.5));
% %       YOUT3 = yout(ind31:ind32,:); TOUT3 = tout(ind31:ind32); TOUT3 = TOUT3-TOUT3(1);
     NADHsim = YOUT2(:,13);
%      [score,scalary] = findSquares(TOUT,NADHsim,dpNADH);
line(TOUT2, NADHsim, 'Color', 'k', 'LineWidth', 2);
%  [score,scalary] = findSquares(TOUT2,YOUT2(:,13),dpNADH);
%      A = find(TOUT2<=2); indA = A(end);
%      [score,scalary] = findSquares(TOUT2(1:indA),NADHsim(1:indA),dpNADH);
%      
%      line(TOUT(1:indA), NADHsim(1:indA), 'Color', 'r', 'LineWidth', 2);
%     hold on;
%      scatter(dpNADH(:,1),dpNADH(:,2).*scalary, 'r','ks');
     hold on
      [score,scalary] = findSquares(TOUT2,YOUT2(:,13),dpNADH2);
      scatter(dpNADH2(:,1),dpNADH2(:,2).*scalary, 100,'r','ko', 'LineWidth', 2);
      h = legend( 'simulation', 'Pan et al. 2019', 'Location', 'North','fontsize',fontsize1);%ref69
     xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('NADH','fontsize',fontsize1);
     box on;
%      
%      
 figure()%%%%%%%%%%%%%%%%%%%%%ATP
     set(gca,'FontName','Times New Roman','FontSize',fontsize2);
     ATPsim = Atot-YOUT2(:,6);
     line(TOUT2, ATPsim, 'Color', 'k', 'LineWidth', 2);
     hold on;
     [score,scalary] = findSquares(TOUT2,ATPsim,dpATP);
     scatter(dpATP(:,1),dpATP(:,2).*scalary,100,'r','o','filled', 'LineWidth', 2);
     hold on;
        [score,scalary] = findSquares(TOUT2,ATPsim,dpATP4);
     scatter(dpATP4(:,1),dpATP4(:,2).*scalary, 100, 'r','o', 'LineWidth', 2);
     xlabel('Time after hypoxia (h)','fontsize',fontsize1);ylabel('ATP','fontsize',fontsize1);
h = legend( 'simulation', 'Casey et al. 2002', 'Wang et al. 1998', 'Location', 'North','fontsize',fontsize1);%ref76;86
% h = legend( 'simulation', 'Casey et al. 2002', 'Location', 'North','fontsize',fontsize1);%ref76;86
box on
  
%
     
% 
     figure()%%%%%%%%%%%%%%%%%%%%%%%%%%%%SIRT1
     set(gca,'FontName','Times New Roman','FontSize',fontsize2);
     SIRT1sim = YOUT2(:,14);
     [score,scalary] = findSquares(TOUT2,SIRT1sim,dpSIRT1_1);
     line(TOUT2, SIRT1sim, 'Color', 'k', 'LineWidth', 2);
     hold on;
     scatter(dpSIRT1_1(:,1),dpSIRT1_1(:,2).*scalary, 100,'r','o','filled');
     hold on;
     [score,scalary] = findSquares(TOUT2,SIRT1sim,dpSIRT3);
     scatter(dpSIRT3(:,1),dpSIRT3(:,2).*scalary, 100,'r','o','LineWidth', 2');
%       scatter(dpSIRT1_2(:,1),dpSIRT1_2(:,2).*scalary, 50,'ko','LineWidth', 2);
      xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('SIRT1','fontsize',fontsize1);
      h = legend( 'simulation', 'Chen et al. 2011', 'Chen et al. 2025', 'Location', 'North','fontsize',fontsize1);%r3-ref37
       box on
% 

% dpHIF1OHP
 figure()%%%%%%%%%%%%%%%%%%%%%%%%%%%%hifoh
     set(gca,'FontName','Times New Roman','FontSize',fontsize2);
     HIFOHsim = YOUT2(:,9);
     [score,scalary] = findSquares(TOUT2,HIFOHsim,dpHIF1OHP);
     line(TOUT2, HIFOHsim, 'Color', 'k', 'LineWidth', 2);
     hold on;
     scatter(dpHIF1OHP(:,1),dpHIF1OHP(:,2).*scalary, 100,'r','o','LineWidth', 2);
%      hold on;
%       scatter(dpSIRT1_2(:,1),dpSIRT1_2(:,2).*scalary, 50,'ko','LineWidth', 2);
      xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('HIFOH','fontsize',fontsize1);
      h = legend( 'simulation', 'Tian et al. 2011',  'Location', 'North','fontsize',fontsize1);%r133
       box on
%       
figure()%%%%%%%%%%%%%%%%%%%%%%%%%%%%HIF
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
HIF1ACsim = YOUT2(:,8);
% HIF1ACsim = YOUT2(:,7)+ YOUT2(:,8)+ YOUT2(:,9);
[score,scalary] = findSquares(TOUT2,HIF1ACsim,dpHIF1a_1);
[score2,scalary2] = findSquares(TOUT2,HIF1ACsim,dpHIF1a2);
% [score3,scalary3] = findSquares(TOUT2,HIF1ACsim,dpHIF1a3);
line(TOUT2, HIF1ACsim, 'Color', 'k', 'LineWidth', 2);
% HIF1atotsim = YOUT2(:,7)+YOUT2(:,8)+YOUT2(:,9);
% [score,scalary] = findSquares(TOUT2,HIF1atotsim,dpHIF1a);
% line(TOUT2, HIF1atotsim, 'Color', 'r', 'LineWidth', 2);

hold on;
scatter(dpHIF1a2(:,1),dpHIF1a2(:,2).*scalary2, 100,'r','o','LineWidth', 2);
% hold on;
% scatter(dpHIF1a3(:,1),dpHIF1a3(:,2).*scalary3, 100,'r','o','filled');
hold on;
scatter(dpHIF1a_1(:,1),dpHIF1a_1(:,2).*scalary, 100,'r','o','filled');
% hold on;
% scatter(dpHIF1a_2(:,1),dpHIF1a_2(:,2).*scalary, 100,'r','o','LineWidth', 2);

xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('HIF1\alphaAC','fontsize',fontsize1);
% xlabel('Time after hypoxia (h)','fontsize',fontsize1); ylabel('total HIF1\alpha','fontsize',fontsize1);
% h = legend( 'simulation','Haugen et al. 2012', 'Cavadas et al. 2015', 'Cavadas et al. 2015', 'Cavadas et al. 2015', 'Location', 'North','fontsize',fontsize1);%ref108;rg11
h = legend( 'simulation','Haugen et al. 2012', 'Cavadas et al. 2015', 'Location', 'North','fontsize',fontsize1);%ref108;rg11
box on




elseif strcmp(Setting,'HighpAMPK1')
    figure()
set(gca,'FontName','Times New Roman','FontSize',16);
    NADratio = YOUT2(:, 12)./YOUT2(:,13); 
    [score,scalary] = findSquares(TOUT2,NADratio,dpNADratio)
    A = find(TOUT2<=2); indA = A(end);
    line(TOUT2, NADratio, 'Color', 'r', 'LineWidth', 2);
    hold on;
     scatter(dpNADratio(:,1),dpNADratio(:,2).*scalary, 'ks','filled');
     xlabel('Time (h)'); ylabel('NAD^{+}/NADH');
     
      figure()
set(gca,'FontName','Times New Roman','FontSize',16);
    pAMPKsim = YOUT2(:,2); pAMPKfold = pAMPKsim./pAMPKsim(1);
    line(TOUT2, pAMPKfold, 'Color', 'r', 'LineWidth', 2);
%      h = legend( 'sim', 'data', 'Location', 'North');
     xlabel('Time (h)'); ylabel('pAMPK fold');
     
     
end
    


% end




