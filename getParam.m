% function [param] = getParam(ver,setting)
function [param] = getParam() %version
 %250605
    k_input_O2 = 0.3; %1;%normal 21% =
     k_etc_f = 3;
     k_etc_b = 0.3; %30
     Jm_etc_O2 = 10;
     k_O2_ROS = 0.5; 
     Jm_O2_ROS = 0.5; 
     k_NOX_f = 5;
     k_NOX_b = 1; %30; %参考R43-REF83 38/s
     Ji_NOX_pAMPK = 1; 
     Jm_NOX_O2 = 10;   
     kd_ROS =10; 
     Ja_SCAV_pAMPK = 0.05;
     ks1_SCAV = 0.1;
     ks2_SCAV = 0.5; 
     Jm_SCAV = 0.2; 
     kd_SCAV = 0.05* kd_ROS; %r3
     k_phos1_AMPK = 0.5; 
     k_phosLKB1_AMPK = 0.5;
     k_phosCaM_AMPK = 0.5; 
     Ja_pAMPK_Aratio = 0.03; %参考r48
     Ja_pAMPK_NADratio = 0.5;
     k_unphos_pAMPK = 1;
     k_a_f = 0.1; 
     k_a_b = 1;
     k_a_f2 = 0.1; 
     k_a_b2 = 1;
     ks1_free = 1;
     ks2_free = 0.5; 
     Ja_HIF1 = 1; 
     kd_free = 0.5; 
     k_AC_free_SIRT1 = 1; 
     k_free_AC_P300 = 0.01;
     k_free_OH_PHD = 0.1; 
     Jm_free_OH_O2 = 0.1;
     k_bind_HIF = 1; 
     k_unbind_HIF = 0.1; 
     kd_OH = 0.5; 
     k_bind_SIRT1 = 0.05; 
     k_unbind_SIRT1 = 1;
     k_NAM_f = 1;
     k_NAM_b = 0.5; 
     Ja_NAM_b_SIRT1 = 0.01;%0.05;
     Ja_NAM_f_pAMPK = 0.05; 
     k_NAD_f = 5;
     k_NAD_b = 0.14;    
     ks1_SIRT1 = 0.02;
     ks2_SIRT1 = 0.2;%%%
     kd_SIRT1 = 0.04;
     KmutantNAD = 0;
      Ja_SIRT1_HIF1a = 0.1;
      Jnad_sirt=10-4;


param = [k_input_O2  k_etc_f  k_etc_b   Jm_etc_O2   k_O2_ROS   Jm_O2_ROS...
   k_NOX_f   k_NOX_b   Ji_NOX_pAMPK   Jm_NOX_O2   kd_ROS   Ja_SCAV_pAMPK   ks1_SCAV...
   ks2_SCAV   Jm_SCAV   kd_SCAV   k_phos1_AMPK   k_phosLKB1_AMPK   k_phosCaM_AMPK...
   Ja_pAMPK_Aratio   Ja_pAMPK_NADratio   k_unphos_pAMPK   k_a_f   k_a_b...
   k_a_f2   k_a_b2   ks1_free   ks2_free   Ja_HIF1   kd_free   k_AC_free_SIRT1...
   Ja_SIRT1_HIF1a   k_free_AC_P300   k_free_OH_PHD   Jm_free_OH_O2   k_bind_HIF...%31
   k_unbind_HIF   kd_OH   k_bind_SIRT1   k_unbind_SIRT1   k_NAM_f   k_NAM_b...
   Ja_NAM_b_SIRT1  Ja_NAM_f_pAMPK   k_NAD_f   k_NAD_b   ks1_SIRT1    ks2_SIRT1...
   kd_SIRT1   KmutantNAD  Jnad_sirt];

