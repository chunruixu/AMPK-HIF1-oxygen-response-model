function dydt = ODE_ROSAMPK(t,y,param)

% ODEs num = 14  %%0903
dydt = zeros(14, 1);
O2 = y(1);
pAMPK = y(2);
ROS = y(3);
SCAV = y(4);%ROS scavenger
deltaH = y(5);
AMP = y(6);
HIF1a_free = y(7);
HIF1a_AC = y(8);
HIF1a_OH = y(9);
HIF1 = y(10);
NAM = y(11);
NAD = y(12);
%R47-in vivo mice bone cell - NAD+ ~15-20pmol/L ~0.02uM
%R13 physiological intracellular NADH 1~10uM
NADH = y(13);
SIRT1 = y(14);

% algebraic - mM


global Atot NAtot AMPKtot
Atot = 30;% 30e3; %[R31] TABLE3; 30mM
NAtot = 2*2;%15; %[R31] TABLE3 3mM;   %[R13] physiological intracellular NADH: 1-10uM=0.001-0.01mM   %NAM in mammalian tissues - R7 - 11~400uM
AMPKtot = 0.1;%3.2e-4mM;%20ppm - pax database, integrated whole organism
%https://www.sinobiological.com/resource/ampk/proteins#:~:text=Molecular%20weight,%28Da%29%2062319.61
%mass: ~62000Da

ATP = Atot - AMP;  %R-M2-REF136 ADP:ATP~0.05
SIRT1__NAM = NAtot - NADH - NAD - NAM;
%R7REF58 11-500uM
AMPK = AMPKtot - pAMPK;

n1 = 1; n2 = 2; n3= 1;
% initial condition

%parameters
k_input_O2 = param(1); k_etc_f = param(2); k_etc_b = param(3);Jm_etc_O2 = param(4);
k_O2_ROS = param (5); Jm_O2_ROS = param(6); k_NOX_f = param(7); k_NOX_b = param(8); 
Ji_NOX_pAMPK = param(9); Jm_NOX_O2 = param(10);   
kd_ROS = param(11); Ja_SCAV_pAMPK = param(12);
ks1_SCAV = param(13); ks2_SCAV = param(14); Jm_SCAV = param(15); kd_SCAV = param(16);
k_phos1_AMPK = param(17); k_phosLKB1_AMPK = param(18); k_phosCaM_AMPK = param(19); 
Ja_pAMPK_ARatio = param(20);Ja_pAMPK_NADRatio = param(21); k_unphos_pAMPK = param(22);
k_a_f = param(23); k_a_b = param(24); k_a_f2 = param(25); k_a_b2 = param(26);
ks1_free = param(27); ks2_free = param(28); Ja_HIF1 = param(29); kd_free = param(30); k_AC_free_SIRT1 = param(31); 
Ja_SIRT1_HIF1a = param(32); k_free_AC_P300 = param(33); k_free_OH_PHD = param(34); Jm_free_OH_O2 = param(35);
k_bind_HIF = param(36); k_unbind_HIF = param(37); kd_OH = param(38); k_bind_SIRT1 = param(39); k_unbind_SIRT1 = param(40);
k_NAM_f = param(41); k_NAM_b = param(42); Ja_NAM_b_SIRT1 = param(43); Ja_NAM_f_pAMPK = param(44); 
k_NAD_f = param(45); k_NAD_b = param(46); ks1_SIRT1 = param(47); ks2_SIRT1 = param(48); kd_SIRT1 = param(49);

% alpha = param(50); 
KmutantNAD = param(50);%for adding NAD - RA10
Jnad_sirt=param(51);

%% 230614 ROS-AMPK equations
ARatio = AMP/ATP;
NADRatio = NAD/NADH;
ETC_F = k_etc_f*NADH*O2^0.5/(Jm_etc_O2^0.5+O2^0.5); %forward reaction 1
ETC_B = k_etc_b*NAD*deltaH^10; %backward reaction 1
% n1 = 1; %reaction 3

O2_ROS = k_O2_ROS*O2^n1/(Jm_O2_ROS^n1+O2^n1);
% n2 = 2;%1; %reaction 4
% NOX_F = k_NOX_f*(Ji_NOX_pAMPK+alpha*pAMPK)/(Ji_NOX_pAMPK+pAMPK)*NADH*O2^n2/(Jm_NOX_O2^n2+O2^n2); %forward reaction 4
NOX_F = k_NOX_f*NADH*O2^n2/(Jm_NOX_O2^n2+O2^n2); %forward reaction 4%%%250608

% NOX_B = 0*k_NOX_b*NAD*ROS;%%250608去掉或极大减弱ROS-nox-NADH这个方向的反应
% AMPK_P = (k_ampk_phos1 + k_ampk_phosLKB*ARatio + k_ampk_phosCaM*ROS/(ROS+Ja_rosampk))*AMPK;%reaction 5
AMPK_P = (k_phos1_AMPK + k_phosLKB1_AMPK*ARatio/(ARatio+Ja_pAMPK_ARatio) + k_phosCaM_AMPK*NADRatio/(NADRatio+Ja_pAMPK_NADRatio))*AMPK;%reaction 5
AMPK_UP = k_unphos_pAMPK*pAMPK;


A_F = k_a_f*AMP*deltaH^3; %reaction 2
A_B = k_a_b*ATP;
A_F2 = k_a_f2*AMP;
A_B2 = k_a_b2*ATP;



AC_FREE_SIRT1 = k_AC_free_SIRT1*SIRT1/(SIRT1+Ja_SIRT1_HIF1a)*HIF1a_AC*NAD;
%0428
% AC_FREE_SIRT1 = k_AC_free_SIRT1*SIRT1*HIF1a_AC*NAD;

FREE_AC_P300 = k_free_AC_P300*HIF1a_free*NAM;
FREE_OH_PHD = k_free_OH_PHD*O2^n3/(O2^n3+Jm_free_OH_O2)*HIF1a_free;

HIF_B = k_bind_HIF*HIF1a_AC;
HIF_UB = k_unbind_HIF*HIF1;
SIRT1_B = k_bind_SIRT1*NAM*SIRT1;
SIRT1_UB = k_unbind_SIRT1*SIRT1__NAM;

NAM_F = k_NAM_f*pAMPK/(pAMPK+Ja_NAM_f_pAMPK)*NAM; %0905
NAM_B = k_NAM_b*SIRT1/(SIRT1+Ja_NAM_b_SIRT1)*NAD;


NAD_F = k_NAD_f*NADH;
NAD_B = k_NAD_b*NAD;

%ODEs

dO2dt = k_input_O2 - 1/2 * ETC_F + 1/2 * ETC_B - O2_ROS - NOX_F - FREE_OH_PHD;

dpAMPKdt = AMPK_P - AMPK_UP;

dROSdt = O2_ROS + NOX_F - kd_ROS * pAMPK/(pAMPK + Ja_SCAV_pAMPK) * SCAV * ROS;

dSCAVdt = ks1_SCAV + ks2_SCAV * ROS^4/(ROS^4 + Jm_SCAV^4) - kd_SCAV * SCAV;

ddeltaHdt = 10*ETC_F - 10*ETC_B - 3*A_F + 3*A_B;

dAMPdt = -A_F + A_B - A_F2 + A_B2;

dHIF1a_freedt = ks1_free + ks2_free*HIF1/(Ja_HIF1+HIF1) - kd_free * HIF1a_free + AC_FREE_SIRT1 - FREE_AC_P300 - FREE_OH_PHD;

dHIF1a_ACdt = - AC_FREE_SIRT1 + FREE_AC_P300 - HIF_B + HIF_UB;

dHIF1a_OHdt = FREE_OH_PHD - kd_OH*HIF1a_OH;

dHIF1dt = HIF_B - HIF_UB;

dNAMdt = - SIRT1_B + SIRT1_UB - NAM_F +NAM_B - FREE_AC_P300 + AC_FREE_SIRT1;

% dNADdt = KmutantNAD + ETC_F - ETC_B + NOX_F + NAM_F - NAM_B + NAD_F - NAD_B + FREE_AC_P300 - AC_FREE_SIRT1;
dNADdt = 0 + ETC_F - ETC_B + NOX_F + NAM_F - NAM_B + NAD_F - NAD_B + FREE_AC_P300 - AC_FREE_SIRT1;

dNADHdt = -ETC_F + ETC_B - NOX_F -NAD_F + NAD_B;

% dSIRT1dt = ks1_SIRT1 +ks2_SIRT1 * HIF1/(HIF1 + Ja_SIRT1_HIF1a) - kd_SIRT1 + SIRT1_UB - SIRT1_B;
%0407
dSIRT1dt = ks1_SIRT1 +ks2_SIRT1 * NAD/(NAD+Jnad_sirt) - kd_SIRT1*SIRT1 + SIRT1_UB - SIRT1_B;%change250604 - r3ref36

dydt(1) = dO2dt; dydt(2) = dpAMPKdt; dydt(3) = dROSdt; dydt(4) = dSCAVdt; dydt(5) = ddeltaHdt;...
dydt(6) = dAMPdt; dydt(7) = dHIF1a_freedt; dydt(8) = dHIF1a_ACdt; dydt(9) = dHIF1a_OHdt;...
dydt(10) = dHIF1dt; dydt(11) = dNAMdt; dydt(12) = dNADdt; dydt(13) = dNADHdt; dydt(14) = dSIRT1dt;

end

