close all
clear all

mutant = 'HighSIRT1';
mutant = 'HighpAMPK1';
if strcmp(mutant,'HighpAMPK1')

%pAMPK
% kp = [0.1ALL, 0.5ALL, 0.1, 0.5, 1, 2, 5, 10, 20, 50];%0.1ALL, 0.5ALL
pAMPK = [1.9258e-04, 2.9048e-04, 3.4961e-04,  3.5546e-04, 3.6329e-04,  3.8078e-04, 4.4797e-04, 5.9663e-04, 9.4470e-04,  0.0020];
pAMPKfoldss = [0.5301, 0.7996, 0.9663, 0.9784, 1, 1.0481, 1.2331, 1.6423,  2.6004,  5.6100];
% pAMPKpeak
HIFACfoldss = [1.1050,1.0323, 1.0110, 1.0072, 1, 0.9823, 0.9054, 0.7475, 0.5152, 0.2562];
HIFACss = [1.4956e-04, 1.3972e-04, 1.3700e-04, 1.3633e-04, 1.3535e-04, 1.3296e-04, 1.2255e-04, 1.0117e-04, 6.9728e-05,  3.4679e-05];
HIF1ss = [ 7.9484e-07, 7.4259e-07, 7.2812e-07, 7.2452e-07, 7.1934e-07, 7.0661e-07, 6.5130e-07, 5.3770e-07, 3.7058e-07, 1.8431e-07];
NADratio = [33.5603, 8.7155, 4.5306, 4.3393, 4.1124, 3.7007, 2.7919, 2.0461, 1.5643, 1.2796];
SIRT1 = [0.3698, 0.4119, 0.4191, 0.4195, 0.4201, 0.4212, 0.4243, 0.4280, 0.4318, 0.4351];
elseif strcmp(mutant,'HighSIRT1')
%SIRT activator
ks = [0.1, 0.5, 1, 2, 5, 10, 20, 50];
SIRT1 = [0.4117, 0.4154,0.4201,0.4294, 0.4573,0.5039,0.5974, 0.8786];
SIRT1fold = [ 0.9801,0.9889,1,1.0221,1.0886,1.1996,1.4220,2.0914];
pAMPK = [3.5904e-04,3.6093e-04,3.6329e-04,3.6804e-04,3.8244e-04, 4.0693e-04,4.5726e-04,6.1332e-04];
pAMPKfoldss = pAMPK./pAMPK(3);
HIFACfoldss = [1.0178,1.0098,1,0.9808, 0.9266,0.8462,0.7162, 0.4786];
HIFACss = [ 1.3776e-04, 1.3668e-04, 1.3535e-04,1.3275e-04,1.2541e-04,1.1454e-04,9.6944e-05,6.4775e-05];
HIF1ss = [7.3212e-07, 7.2640e-07,7.1934e-07, 7.0554e-07,6.6653e-07,6.0874e-07,5.1522e-07,3.4426e-07];
NADratio = [4.0261,4.0644,4.1124,4.2093, 4.5050,5.0134,6.0765, 9.4589];

[R3, P3] = corrcoef(HIFACss, SIRT1fold);
end


%% R and P
[R1, P1] = corrcoef(pAMPKfoldss, HIFACss);
[R2, P2] = corrcoef(pAMPK, HIFACss);


%% PLOTS
fontsize1 = 22; fontsize2 = 18;

figure()%HIFAC VS pAMPK
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
scatter(pAMPKfoldss,HIFACss,50,'k','o','filled');
hold on;
line(pAMPKfoldss, HIFACss, 'Color', 'r','Linestyle', '--', 'LineWidth', 2);
xlabel('pAMPK Change Fold','fontsize',fontsize1); ylabel('HIF1\alpha-AC','fontsize',fontsize1);
% ylim([0.8*min(HIFACss) 1.25*max(HIFACss)]);  xlim([0.8*min(pAMPKfoldss) 1.25*max(pAMPKfoldss)]);
box on 

figure()%SIRT1 VS pAMPK
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
scatter(pAMPKfoldss,SIRT1,50,'k','o','filled');
hold on;
line(pAMPKfoldss, SIRT1, 'Color', 'r','Linestyle', '--', 'LineWidth', 2);
xlabel('pAMPK Change Fold','fontsize',fontsize1); ylabel('SIRT1','fontsize',fontsize1);
box on 

figure()%nadratio VS pAMPK
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
scatter(pAMPKfoldss,NADratio,50,'k','o','filled');
hold on;
line(pAMPKfoldss, NADratio, 'Color', 'r','Linestyle', '--', 'LineWidth', 2);
xlabel('pAMPK Change Fold','fontsize',fontsize1); ylabel('NAD^{+}/NADH','fontsize',fontsize1);
box on 

if strcmp(mutant,'HighSIRT1')
figure()%HIF VS SIRT1
set(gca,'FontName','Times New Roman','FontSize',fontsize2);
scatter(SIRT1fold,HIFACss,50,'k','o','filled');
hold on;
line(SIRT1fold,HIFACss, 'Color', 'r','Linestyle', '--', 'LineWidth', 2);
num_ticks = 5; % 指定要显示的刻度标签数目
xtick_pos = linspace(min(SIRT1fold), max(SIRT1fold), num_ticks); % 生成等间距的刻度位置
% 设置 x 轴的刻度位置和标签
xticks(xtick_pos);
% xticklabels(cellstr(num2str(xtick_pos', '%.0f'))); % 将刻度位置转换为字符串标签

xlabel('SIRT1 Change Fold','fontsize',fontsize1); ylabel('HIF1\alpha-AC','fontsize',fontsize1);
box on 
end
% figure()%HIF1 VS pAMPK
% set(gca,'FontName','Times New Roman','FontSize',fontsize2);
% scatter(pAMPKfoldss,HIF1ss,50,'k','o','filled');
% hold on;
% line(pAMPKfoldss, HIF1ss, 'Color', 'r','Linestyle', '--', 'LineWidth', 2);
% xlabel('pAMPK Change Fold','fontsize',fontsize1); ylabel('HIF1','fontsize',fontsize1);
% box on 
% 
% figure()%HIFAC VS pAMPK
% set(gca,'FontName','Times New Roman','FontSize',fontsize2);
% scatter(pAMPK,HIFACss,50,'k','o','filled');
% hold on;
% line(pAMPK, HIFACss, 'Color', 'r','Linestyle', '--', 'LineWidth', 2);
% xlabel('pAMPK','fontsize',fontsize1); ylabel('HIF1\alpha-AC','fontsize',fontsize1);
% box on 