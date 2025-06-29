function [param] = getMutantParam(param,setting)
%根据给定的参数和设置（setting）生成突变后的参数

%从输入的原始参数 param 中提取各个参数的数值
k_input_O2 = param(1); k_etc_f = param(2); k_etc_b = param(3);Jm_etc_O2 = param(4);
k_O2_ROS = param(5); Jm_O2_ROS = param(6); k_NOX_f = param(7); k_NOX_b = param(8); 
Ji_NOX_pAMPK = param(9); Jm_NOX_O2 = param(10);   
kd_ROS = param(11); Ja_SCAV_pAMPK = param(12);
ks1_SCAV = param(13); ks2_SCAV = param(14); Jm_SCAV = param(15); kd_SCAV = param(16);
k_phos1_AMPK = param(17); k_phosLKB1_AMPK = param(18); k_phosCaM_AMPK = param(19); 
Ja_pAMPK_Aratio = param(20);Ja_pAMPK_NADratio = param(21); k_unphos_pAMPK = param(22);
k_a_f = param(23); k_a_b = param(24); k_a_f2 = param(25); k_a_b2 = param(26);
ks1_free = param(27); ks2_free = param(28); Ja_HIF1 = param(29); kd_free = param(30); k_AC_free_SIRT1 = param(31); 
Ja_SIRT1_HIF1a = param(32); k_free_AC_P300 = param(33); k_free_OH_PHD = param(34); Jm_free_OH_O2 = param(35);
k_bind_HIF = param(36); k_unbind_HIF = param(37); kd_OH = param(38); k_bind_SIRT1 = param(39); k_unbind_SIRT1 = param(40);
k_NAM_f = param(41); k_NAM_b = param(42); Ja_NAM_b_SIRT1 = param(43); Ja_NAM_f_pAMPK = param(44); 
k_NAD_f = param(45); k_NAD_b = param(46); ks1_SIRT1 = param(47); ks2_SIRT1 = param(48); kd_SIRT1 = param(49);

KmutantNAD = param(50);  %NAtot = param(34);
Jnad_sirt=param(51);

%两种设置：'hypoxia1' 和 'HighpAMPK1'
    if strcmp(setting,'hypoxia1')
        %比较字符串,见MD.2。如果二者相同，则返回 1 (true)，否则返回 0
%         fold = 1/7;%10%   1% 5% 10% 20% 50%
%           fold = 0.05;%20%
          fold = 1/5;
        k_input_O2 = k_input_O2*fold;%氧气输入速率 k_o2input 减小为原来的 1/10
        elseif strcmp(setting,'hypoxia0')
          fold = 0.5;
        k_input_O2 = k_input_O2*fold;%氧气输入速率 k_o2input 减小为原来的 1/2
    elseif strcmp(setting,'hypoxia2')
          fold = 0.1;
        k_input_O2 = k_input_O2*fold;%氧气输入速率 k_o2input 减小为原来的 1/10
           elseif strcmp(setting,'hypoxia3')
          fold = 0.05;
        k_input_O2 = k_input_O2*fold;
        elseif strcmp(setting,'hyperoxia0')
        fold = 1.1;
%         fold = 5;
        k_input_O2 = k_input_O2*fold;
    elseif strcmp(setting,'hyperoxia1')
        fold = 1.5;
%         fold = 5;
        k_input_O2 = k_input_O2*fold;
         elseif strcmp(setting,'hyperoxia2')
        fold = 2;
        k_input_O2 = k_input_O2*fold;%
        elseif strcmp(setting,'hyperoxia3')
        fold = 5;
        k_input_O2 = k_input_O2*fold;%
    elseif strcmp(setting,'HighpAMPK1')
        fold = 2;
        k_phos1_AMPK = k_phos1_AMPK*fold;
%         k_phosLKB1_AMPK = k_phosLKB1_AMPK*fold;
%         k_phosCaM_AMPK  = k_phosCaM_AMPK*fold;
        %将 AMPK 的磷酸化速率相关的三个参数 k_ampk_phos1、k_ampk_phosLKB 和 k_ampk_phosCaM 增大为原来的 2 倍
           elseif strcmp(setting,'HighpAMPKALL')
        fold = 0.1;
        k_phos1_AMPK = k_phos1_AMPK*fold;
        k_phosLKB1_AMPK = k_phosLKB1_AMPK*fold;
        k_phosCaM_AMPK  = k_phosCaM_AMPK*fold;
    elseif strcmp(setting,'AMPKinhibitor&hypoxia')
 k_input_O2 = k_input_O2/7;%氧气输入速率 k_o2input 
  k_phos1_AMPK = k_phos1_AMPK*0.2;
    elseif strcmp(setting,'ADPATP')
        fold = 2;
        k_a_f2 =fold*k_a_f2;
         elseif strcmp(setting,'NMN')
        KmutantNAD = 2;
  elseif strcmp(setting,'HighSIRT1')
        fold = 2;
        ks1_SIRT1 =fold*ks1_SIRT1;
    end

    %据设置后的参数，更新参数向量 param
    param = [k_input_O2;k_etc_f; k_etc_b; Jm_etc_O2;k_O2_ROS; Jm_O2_ROS;...
            k_NOX_f;k_NOX_b;Ji_NOX_pAMPK; Jm_NOX_O2;   kd_ROS; Ja_SCAV_pAMPK;ks1_SCAV;...
            ks2_SCAV; Jm_SCAV; kd_SCAV;k_phos1_AMPK; k_phosLKB1_AMPK;k_phosCaM_AMPK;...
            Ja_pAMPK_Aratio; Ja_pAMPK_NADratio;k_unphos_pAMPK;k_a_f; k_a_b;...
            k_a_f2; k_a_b2;ks1_free;ks2_free; Ja_HIF1; kd_free; k_AC_free_SIRT1;...
            Ja_SIRT1_HIF1a;k_free_AC_P300;k_free_OH_PHD; Jm_free_OH_O2;k_bind_HIF;...
            k_unbind_HIF; kd_OH; k_bind_SIRT1; k_unbind_SIRT1;k_NAM_f;k_NAM_b;...
            Ja_NAM_b_SIRT1;Ja_NAM_f_pAMPK; k_NAD_f;k_NAD_b;ks1_SIRT1;ks2_SIRT1;...
              kd_SIRT1;KmutantNAD; Jnad_sirt];
end