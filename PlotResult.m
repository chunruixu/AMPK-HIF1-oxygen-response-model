function [] = PlotResult(tout, yout, Time1, Time2)

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
% 这块不会
SIRT1__NAM = NAtotyout - yout(:, 11) - yout(:, 12) - yout(:, 13);
Atotyout = zeros(timeL, 1) + Atot;
ATP = Atotyout - yout(:,6);

%% time course plot %%


ind1 = find(tout>=Time1); ind2 = find(tout>=Time2);
TOUT = tout(ind1:ind2); YOUT = yout(ind1:ind2,:);

figure()%O2, ROS
 p1 = line(TOUT, YOUT(:, 1), 'Color', 'r', 'LineWidth', 2, 'Linestyle', ':');
  p3 = line(TOUT, YOUT(:, 6), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-')
  p4 = line(TOUT, ATP(ind1:ind2), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '--')
  h = legend('O2', 'ADP','ATP', 'Location', 'North','fontsize',fontsize1);
hold on
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
     
     
     figure()%O2, ROS
 p2 = line(TOUT, YOUT(:, 3), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
  h = legend( 'ROS', 'Location', 'North','fontsize',fontsize1);
hold on
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);

  figure()%NAD & NADH & NAM
 p5 = line(TOUT, YOUT(:, 13), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
 p6 = line(TOUT, YOUT(:, 12), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '--');
 p7 = line(TOUT, YOUT(:, 11), 'Color', 'r', 'LineWidth', 2, 'Linestyle', ':');
%  p8 = line(TOUT, SIRT1__NAM(ind1:ind2), 'Color', 'g', 'LineWidth', 2, 'Linestyle', '-');
  p7 = line(TOUT, YOUT(:, 14), 'Color', 'g', 'LineWidth', 2, 'Linestyle', ':');
hold on
% xa = [Insert 3.5];
% ya = [Insert 3];
% drawarrow(xa,ya,'textarrow',gca,'Hypoxia','b');%drawarrow(startpoint,endpoint,'String','color',linewidth);

%  h = legend( 'NADH', 'NAD+', 'NAM','deltaH', 'Location', 'North');
  h = legend( 'NADH', 'NAD^+','NAM', 'SIRTfree', 'Location', 'North','fontsize',fontsize1);
%   h = legend( 'NADH', 'NAD+','deltaH', 'Location', 'North');
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);


 figure()%AMPK & ratio
  p9 = line(TOUT, YOUT(:, 2), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
   p11 = line(TOUT, YOUT(:, 12)./YOUT(:, 13), 'Color', 'b', 'LineWidth', 2, 'Linestyle', '--');
  h = legend( 'pAMPK', 'NAD+/NADH', 'Location', 'North','fontsize',fontsize1);
  
   figure()%AMPK
  p9 = line(TOUT, YOUT(:, 2), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
  h = legend( 'pAMPK',  'Location', 'North','fontsize',fontsize1);

 figure()%RATIO1
 
  p12 = line(TOUT, YOUT(:,6)./ATP(ind1:ind2), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
  h = legend( 'AMP/ATP', 'Location', 'North','fontsize',fontsize1);
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);


 figure()%RATIO2

   p13 = line(TOUT, YOUT(:, 12)./YOUT(:, 13), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
  h = legend( 'NAD+/NADH', 'Location', 'North','fontsize',fontsize1);
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);
  
 figure()%HIF1a_free ;% HIF1a_AC ;% HIF1a_OH;% HIF1 ;
%   p14 = line(TOUT, YOUT(:, 7), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
%   p15 = line(TOUT, YOUT(:, 8), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '--');
%   p16 = line(TOUT, YOUT(:, 9), 'Color', 'g', 'LineWidth', 2, 'Linestyle', '--');
  p17 = line(TOUT, YOUT(:, 10), 'Color', 'b', 'LineWidth', 2, 'Linestyle', ':');
  h = legend( 'HIF1', 'Location', 'North','fontsize',fontsize1);
%     h = legend( 'HIF1a_free', 'HIF1a_AC',  'HIF1', 'Location', 'North','fontsize',fontsize1);

%   h = legend( 'HIF1a_free', 'HIF1a_AC', 'HIF1a_OH', 'HIF1', 'Location', 'North','fontsize',fontsize1);
%     xa = [Insert 750];
%     ya = [Insert 650];
% drawarrow(xa,ya,'textarrow',gca,'Hypoxia','b');%drawarrow(startpoint,endpoint,'String','color',linewidth);
     xlabel('Time (h)','fontsize',fontsize1);% ylabel('pAMPK fold','fontsize',fontsize1);

  figure()%SIRT1
  p18 = line(TOUT, YOUT(:, 14)+SIRT1__NAM(ind1:ind2), 'Color', 'r', 'LineWidth', 2, 'Linestyle', '-');
   h = legend( 'SIRT1tot', 'Location', 'North','fontsize',fontsize1);
  xlabel('Time (h)','fontsize',fontsize1);% ylabel(
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
   dpNADratio = [0 2.88; 2 3.68; 4 5.29; 8 5.18]; %ref6-ref45, part, AMPK activator
%    dpNAD1 = [0 0.315; 0.17 0.331; 1 0.442; 2 0.496];%R70
   dpNAD = [0 1; 8 0.8; 16 0.5];%R106REF34
dpNADH = [0.083333333 0.52; 0.166666667 0.97; 0.25 1.81; 0.333333333 1.84; 0.416666667 1.82];%R69-参考意义稍微有点小，因为常氧咋也变化啊
dpATP =[0 16.63; 0.5 14.54; 0.75 12.68; 1 10.81];%r76 
dpATP2 = [0 24.77; 0.5 17; 1 11.61; 1.5 1.1];%R78  
dpATP3 = [0 17.99;0.25	13.49;0.5	4.63;0.75	3.41;1	2.41];%R82
dpATP4 = [0	25.76;1.5	18.28;3	14.27;5	8.05];%R86
dpHIF1a_AC = [0	 0.06; 4  1.85; 8  2.3];%ref108: 綜合幾個fig
dpSIRT1 = [0 0.520; 2 0.754; 4 1.254];







  




