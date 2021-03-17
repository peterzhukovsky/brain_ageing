%% DATASETS
% ATX
datanames={'1001_22053';'1002_22515';'1005_22547';'1009_11509';'1013_22570';'1014_22575';'1015_22590';'1016_22589';'1017_22599';'1018_22710';'1020_22740';'1021_22733';'1022_23418';'1023_22768';'1024_22765';'1025_22795';'1026_22831';'1027_22972';'1029_19429';'1030_23307';'1031_23611';'1033_23323';'1034_23538';'1035_23469';'1037_23754';'1038_23620';'1039_23446';'1041_20445';'1042_18737';'1045_23844';'1046_23768';'1047_23784';'1048_23802';'1049_24299';'1050_23823';'1051_24326';'1053_23932';'1054_24005';'1056_24332';'2008_22644';'2010_22656';'2014_22655';'2015_22667';'2016_22683';'2018_22724';'2022_22893';'2027_20647';'2028_23041';'2031_23105';'2034_23160';'2037_23359';'2038_23383';'2040_23346';'2044_23581';'2046_23357';'2050_23415';'2051_23422';'2052_18801';'2053_23414';'2054_23445';'2055_23444';'2056_23591';'2057_23468';'2059_23532';'2060_23545';'2061_24331';'2063_23667';'2065_24296';'2067_23703';'2068_18726';'2071_20615';'2072_23765';'2074_23786';'2075_23793';'2076_23873';'2077_23848';'2078_23853';'2079_23915';'2080_23920'};
% GSK
datanames={'1001_11280';'1002_11367';'1003_11945';'1004_11449';'1005_11469';'1006_11455';'1007_11464';'1008_11865';'1009_12507';'1010_11940';'1011_12938';'1012_13006';'1013_12653';'1014_12367';'1015_13045';'1016_13057';'1017_13070';'1018_13128';'1022_13125';'1023_13289';'1024_13599';'1029_13505';'1050_13560';'1051_12366';'1053_11463';'1054_12506';'1055_13246';'1056_13260';'1057_13280';'2001_10684';'2002_10709';'2003_12306';'2004_11466';'2005_11665';'2006_10801';'2008_11371';'2009_11861';'2010_12481';'2011_12940';'2012_11637';'2013_11973';'2014_11967';'2015_11544';'2016_13200';'2017_11965';'2018_13047';'2020_13274';'2023_13677';'2026_13467';'2029_13356';'2038_13787';'2051_10792';'2052_10704';'2053_10802';'2054_10793';'2055_11260';'2056_11274';'2058_11268';'2059_11269';'2060_11277';'2061_11370';'2062_11535';'2063_11543';'2064_11673';'2065_11661';'2066_11674';'2067_11862';'2068_11913';'2069_11914';'2070_11960';'2071_11928';'2073_12397';'2074_12520';'2075_12663';'2076_12746';'2077_13049';'2078_13056';'2079_13107';'2080_13130';'2081_13143';'2082_13178';'2083_13204';'2084_13205';'2085_13233';'2086_13305';'2087_13355';'2088_13376';'2089_13466';'2090_13469';'2091_13501';'2094_13696';'2095_13697'};
% Ageing
datanames={'1001_20267';'1002_20282';'1003_20284';'1004_11903';'1005_20319';'1006_20317';'1007_20336';'1008_20402';'1009_20399';'1010_20372';'1011_20480';'1012_20445';'1013_20484';'1014_20483';'1015_20504';'1016_20506';'1017_20565';'1018_18209';'1019_20588';'1020_14343';'1021_20586';'1022_20660';'1025_19429';'1026_20754';'1027_20825';'1028_20901';'1029_20900';'1030_20924';'1031_20961';'1032_21006';'1033_21022';'2001_20283';'2002_20325';'2003_20338';'2004_20342';'2005_20397';'2006_20357';'2007_20405';'2008_20476';'2010_20599';'2011_20628';'2012_20507';'2013_20615';'2014_20646';'2015_20647';'2016_20659';'2017_20689';'2018_20703';'2019_20705';'2020_20690';'2021_20726';'2022_20753';'2023_20814';'2025_20813';'2026_20820';'2027_20839';'2028_20868';'2029_20918';'2030_20917';'2031_20931';'2032_20935';'2033_20939';'2034_20957';'2035_20964';'2036_20990'}
%ENDO 
datanames={'1001_11655';'1002_10846';'1003_14296';'1004_16928';'1005_17287';'1006_17291';'1007_14061';'1008_17485';'1009_17577';'1010_17576';'1011_16073';'1012_17636';'1013_17645';'1014_14693';'1015_17694';'1016_17713';'1017_17718';'1018_10058';'1020_17756';'1021_17765';'1022_17777';'1023_17782';'1024_17810';'1025_17828';'1026_17851';'1027_17863';'1028_17891';'1029_17976';'1031_17970';'1032_18093';'1033_18300';'1034_18299';'1035_18307';'1036_18308';'1037_18561';'1038_60';'1039_18683';'1040_18697';'1041_18717';'1042_16561';'1043_18851';'1044_18905';'1045_18899';'1046_18981';'1047_18992';'1049_19086';'1050_19088';'1051_19132';'1052_19143';'2001_16758';'2002_13502';'2003_16821';'2004_11623';'2005_17034';'2006_17290';'2007_17286';'2008_16846';'2009_17295';'2010_17605';'2013_17711';'2014_17708';'2015_17733';'2016_17750';'2017_17861';'2018_17912';'2019_17928';'2020_17920';'2021_17910';'2022_17990';'2024_18070';'2025_18087';'2026_18112';'2027_18080';'2028_18164';'2029_18135';'2030_18199';'2031_18328';'2032_18413';'2033_18446';'2035_18443';'2036_18479';'2037_18528';'2038_18486';'2039_18544';'2040_18548';'2041_18657';'2042_18680';'2043_18705';'2044_18726';'2046_18801';'2047_18838';'2048_18959';'2049_18975';'2050_18998';'2051_19059'};
%ROCK 
datanames1={'A00008326';'A00008399';'A00010893';'A00013809';'A00021039';'A00023510';'A00025566';'A00027544';'A00027651';'A00028150';'A00028152';'A00028177';'A00028184';'A00028185';'A00028192';'A00028207';'A00028246';'A00028266';'A00028287';'A00028339';'A00028340';'A00028352';'A00028380';'A00028389';'A00028399';'A00028400';'A00028429';'A00028430';'A00028468';'A00028552';'A00028605';'A00028606';'A00028613';'A00028625';'A00028656';'A00028678';'A00028691';'A00028694';'A00028753';'A00028754';'A00028766';'A00028784';'A00028822';'A00028842';'A00028844';'A00028845';'A00028912';'A00028929';'A00028994';'A00028995';'A00029075';'A00029076';'A00029092';'A00029104';'A00029126';'A00029127';'A00029215';'A00029219';'A00029230';'A00029231';'A00029303';'A00029304';'A00029979';'A00030912';'A00030947';'A00030989';'A00030990';'A00031166';'A00031167';'A00031272';'A00031452';'A00031459';'A00031549';'A00031604';'A00031605';'A00031625';'A00031683';'A00031794';'A00031871';'A00031872';'A00031893';'A00031894';'A00032007';'A00032010';'A00032862';'A00032875';'A00032876';'A00033011';'A00033021';'A00033232';'A00033247';'A00033248';'A00033589';'A00033609';'A00033640';'A00033643';'A00033647';'A00033673';'A00033734';'A00033735';'A00033747';'A00033748';'A00033752';'A00033764';'A00033832';'A00033834';'A00033849';'A00033903';'A00033963';'A00033998';'A00034093';'A00034193';'A00034385';'A00034400';'A00034827';'A00034854';'A00034987';'A00035036';'A00035072';'A00035363';'A00035378';'A00035437';'A00035440';'A00035504';'A00035535';'A00035553';'A00035562';'A00035606';'A00035608';'A00035625';'A00035633';'A00035699';'A00035744';'A00035765';'A00035798';'A00035827';'A00035840';'A00035879';'A00035881';'A00035926';'A00035943';'A00035951';'A00035960';'A00037110';'A00037111';'A00037112';'A00037266';'A00037267';'A00037377';'A00037378';'A00037386';'A00037396';'A00037421';'A00037439';'A00037445';'A00037476';'A00037492';'A00037510';'A00037511';'A00037522';'A00037531';'A00037545';'A00037582';'A00037588';'A00037635';'A00037646';'A00037647';'A00037725';'A00037726';'A00037767';'A00037768';'A00037783';'A00037784';'A00037817'};
datanames2={'A00037818';'A00037831';'A00037848';'A00037877';'A00038115';'A00038189';'A00038285';'A00038411';'A00038414';'A00038424';'A00038519';'A00038520';'A00038522';'A00038605';'A00038616';'A00038623';'A00038642';'A00038718';'A00038731';'A00038805';'A00038806';'A00038817';'A00038831';'A00038832';'A00038957';'A00038959';'A00038998';'A00039041';'A00039074';'A00039083';'A00039118';'A00039143';'A00039168';'A00039218';'A00039249';'A00039257';'A00039276';'A00039277';'A00039331';'A00039353';'A00039391';'A00039393';'A00039431';'A00039432';'A00039461';'A00039463';'A00039488';'A00039490';'A00039500';'A00039559';'A00039560';'A00039593';'A00039607';'A00039635';'A00039636';'A00039639';'A00039640';'A00039655';'A00039669';'A00039685';'A00039686';'A00039699';'A00039700';'A00039755';'A00039758';'A00039783';'A00039819';'A00039820';'A00039845';'A00039846';'A00039895';'A00039916';'A00039923';'A00039951';'A00039952';'A00039966';'A00039974';'A00040116';'A00040151';'A00040152';'A00040174';'A00040182';'A00040193';'A00040248';'A00040257';'A00040285';'A00040286';'A00040301';'A00040311';'A00040312';'A00040324';'A00040343';'A00040351';'A00040382';'A00040383';'A00040439';'A00040461';'A00040462';'A00040493';'A00040494';'A00040500';'A00040517';'A00040524';'A00040525';'A00040567';'A00040573';'A00040580';'A00040594';'A00040599';'A00040623';'A00040628';'A00040629';'A00040640';'A00040678';'A00040740';'A00040741';'A00040757';'A00040784';'A00040806';'A00040891';'A00040915';'A00040944';'A00041540';'A00043283';'A00043299';'A00043324';'A00043325';'A00043384';'A00043417';'A00043450';'A00043462';'A00043466';'A00043509';'A00043520';'A00043521';'A00043617';'A00043634';'A00043635';'A00043649';'A00043677';'A00043704';'A00043721';'A00043722';'A00043739';'A00043762';'A00043790';'A00043998';'A00044011';'A00044012';'A00044068';'A00044130';'A00044131';'A00044154';'A00044291';'A00044307';'A00044344';'A00044345';'A00044369';'A00044370';'A00044405';'A00044410';'A00044427';'A00045554';'A00045589';'A00045590';'A00050720';'A00050721';'A00050742';'A00050795';'A00050940';'A00050977';'A00050998';'A00051064';'A00051456'};
datanames3={'A00051457';'A00051477';'A00051513';'A00051514';'A00051517';'A00051528';'A00051529';'A00051539';'A00051548';'A00051604';'A00051676';'A00051678';'A00051679';'A00051774';'A00051796';'A00051835';'A00051882';'A00051925';'A00051927';'A00052070';'A00052117';'A00052118';'A00052125';'A00052126';'A00052180';'A00052181';'A00052214';'A00052234';'A00052235';'A00052307';'A00052319';'A00052340';'A00052499';'A00052500';'A00052560';'A00052577';'A00052612';'A00052639';'A00053202';'A00053369';'A00053455';'A00053473';'A00053474';'A00053475';'A00053576';'A00053577';'A00053578';'A00053625';'A00053626';'A00053627';'A00053850';'A00053851';'A00053854';'A00053874';'A00053901';'A00053902';'A00053927';'A00053949';'A00054019';'A00054038';'A00054153';'A00054173';'A00054358';'A00054441';'A00054482';'A00054504';'A00054532';'A00054533';'A00054534';'A00054621';'A00054857';'A00054895';'A00054897';'A00054914';'A00054929';'A00055061';'A00055121';'A00055215';'A00055352';'A00055353';'A00055373';'A00055446';'A00055447';'A00055542';'A00055692';'A00055738';'A00055763';'A00055806';'A00055962';'A00056097';'A00056098';'A00056164';'A00056306';'A00056372';'A00056452';'A00056489';'A00056556';'A00056627';'A00056898';'A00056919';'A00056949';'A00057005';'A00057035';'A00057182';'A00057183';'A00057203';'A00057235';'A00057372';'A00057406';'A00057444';'A00057445';'A00057786';'A00057808';'A00057965';'A00058214';'A00058218';'A00058318';'A00058436';'A00058503';'A00058515';'A00058552';'A00058621';'A00058644';'A00058667';'A00058952';'A00058998';'A00058999';'A00059074';'A00059107';'A00059344';'A00059346';'A00059361';'A00059428';'A00059527';'A00059662';'A00059663';'A00059664';'A00059756';'A00059845';'A00059911';'A00060005';'A00060006';'A00060093';'A00060169';'A00060170';'A00060184';'A00060185';'A00060252';'A00060259';'A00060264';'A00060279';'A00060301';'A00060302';'A00060306';'A00060372';'A00060383';'A00060407';'A00060430';'A00060431';'A00060471';'A00060480';'A00060516';'A00060574';'A00060575';'A00060582';'A00060602';'A00060629';'A00060630';'A00060632';'A00060662';'A00060773';'A00060848';'A00060923';'A00060925'};
datanames4={'A00061001';'A00061021';'A00061203';'A00061204';'A00061210';'A00061263';'A00061274';'A00061275';'A00061276';'A00061281';'A00061284';'A00061387';'A00061483';'A00061485';'A00061597';'A00061598';'A00061600';'A00061634';'A00061647';'A00061656';'A00061666';'A00061709';'A00061711';'A00061725';'A00061727';'A00061728';'A00061767';'A00061787';'A00061790';'A00061806';'A00061881';'A00061882';'A00061962';'A00062210';'A00062248';'A00062255';'A00062266';'A00062268';'A00062282';'A00062285';'A00062288';'A00062292';'A00062351';'A00062368';'A00062411';'A00062416';'A00062917';'A00062920';'A00062923';'A00062926';'A00062928';'A00062929';'A00062934';'A00062937';'A00062942';'A00063000';'A00063008';'A00063031';'A00063103';'A00063110';'A00063326';'A00063368';'A00063424';'A00063450';'A00063455';'A00063470';'A00063473';'A00063479';'A00063481';'A00063589';'A00064053';'A00064061';'A00064081';'A00064323';'A00064328';'A00064415';'A00064580';'A00065302';'A00065379';'A00065478';'A00065480';'A00065487';'A00065572';'A00065574';'A00065617';'A00065717';'A00065722';'A00065743';'A00065749';'A00065790';'A00065873';'A00065935';'A00065974';'A00065991';'A00065995';'A00066013';'A00066087';'A00066091';'A00066130';'A00066153';'A00066154';'A00066162';'A00066232';'A00066236';'A00066237';'A00066245';'A00066246';'A00066282';'A00066302';'A00066319';'A00066389';'A00066419';'A00066460';'A00066539';'A00066735';'A00066781';'A00066788';'A00066800';'A00066812';'A00066820';'A00066822';'A00066827';'A00066831';'A00066864';'A00066865';'A00066871';'A00066926';'A00067174';'A00072203';'A00073230';'A00073283';'A00073296';'A00073308';'A00073330';'A00073525';'A00073529';'A00073600';'A00073611';'A00073677';'A00073705';'A00073942';'A00073953';'A00074000'}

%datanames=vertcat(datanames1, datanames2, datanames3, datanames4)
for subject=1:length(datanames)
    subject_id=datanames(subject); subject_id=subject_id{:};subject_id
% Left side
clear SurfArea GrayVol ThickAvg MeanCurv GausCurv FoldInd CurvInd
llabel_names={'L_10d_ROI.label';'L_10pp_ROI.label';'L_10r_ROI.label';'L_10v_ROI.label';'L_11l_ROI.label';'L_13l_ROI.label';'L_1_ROI.label';'L_23c_ROI.label';'L_23d_ROI.label';'L_24dd_ROI.label';'L_24dv_ROI.label';'L_25_ROI.label';'L_2_ROI.label';'L_31a_ROI.label';'L_31pd_ROI.label';'L_31pv_ROI.label';'L_33pr_ROI.label';'L_3a_ROI.label';'L_3b_ROI.label';'L_43_ROI.label';'L_44_ROI.label';'L_45_ROI.label';'L_46_ROI.label';'L_47l_ROI.label';'L_47m_ROI.label';'L_47s_ROI.label';'L_4_ROI.label';'L_52_ROI.label';'L_55b_ROI.label';'L_5L_ROI.label';'L_5m_ROI.label';'L_5mv_ROI.label';'L_6a_ROI.label';'L_6d_ROI.label';'L_6ma_ROI.label';'L_6mp_ROI.label';'L_6r_ROI.label';'L_6v_ROI.label';'L_7AL_ROI.label';'L_7Am_ROI.label';'L_7m_ROI.label';'L_7PC_ROI.label';'L_7PL_ROI.label';'L_7Pm_ROI.label';'L_8Ad_ROI.label';'L_8Av_ROI.label';'L_8BL_ROI.label';'L_8BM_ROI.label';'L_8C_ROI.label';'L_9-46d_ROI.label';'L_9a_ROI.label';'L_9m_ROI.label';'L_9p_ROI.label';'L_a10p_ROI.label';'L_A1_ROI.label';'L_a24pr_ROI.label';'L_a24_ROI.label';'L_a32pr_ROI.label';'L_a47r_ROI.label';'L_A4_ROI.label';'L_A5_ROI.label';'L_a9-46v_ROI.label';'L_AAIC_ROI.label';'L_AIP_ROI.label';'L_AVI_ROI.label';'L_d23ab_ROI.label';'L_d32_ROI.label';'L_DVT_ROI.label';'L_EC_ROI.label';'L_FEF_ROI.label';'L_FFC_ROI.label';'L_FOP1_ROI.label';'L_FOP2_ROI.label';'L_FOP3_ROI.label';'L_FOP4_ROI.label';'L_FOP5_ROI.label';'L_FST_ROI.label';'L_H_ROI.label';'L_i6-8_ROI.label';'L_IFJa_ROI.label';'L_IFJp_ROI.label';'L_IFSa_ROI.label';'L_IFSp_ROI.label';'L_Ig_ROI.label';'L_IP0_ROI.label';'L_IP1_ROI.label';'L_IP2_ROI.label';'L_IPS1_ROI.label';'L_LBelt_ROI.label';'L_LIPd_ROI.label';'L_LIPv_ROI.label';'L_LO1_ROI.label';'L_LO2_ROI.label';'L_LO3_ROI.label';'L_MBelt_ROI.label';'L_MIP_ROI.label';'L_MI_ROI.label';'L_MST_ROI.label';'L_MT_ROI.label';'L_OFC_ROI.label';'L_OP1_ROI.label';'L_OP2-3_ROI.label';'L_OP4_ROI.label';'L_p10p_ROI.label';'L_p24pr_ROI.label';'L_p24_ROI.label';'L_p32pr_ROI.label';'L_p32_ROI.label';'L_p47r_ROI.label';'L_p9-46v_ROI.label';'L_PBelt_ROI.label';'L_PCV_ROI.label';'L_PeEc_ROI.label';'L_PEF_ROI.label';'L_PFcm_ROI.label';'L_PFm_ROI.label';'L_PFop_ROI.label';'L_PF_ROI.label';'L_PFt_ROI.label';'L_PGi_ROI.label';'L_PGp_ROI.label';'L_PGs_ROI.label';'L_PHA1_ROI.label';'L_PHA2_ROI.label';'L_PHA3_ROI.label';'L_PH_ROI.label';'L_PHT_ROI.label';'L_PI_ROI.label';'L_Pir_ROI.label';'L_PIT_ROI.label';'L_pOFC_ROI.label';'L_PoI1_ROI.label';'L_PoI2_ROI.label';'L_POS1_ROI.label';'L_POS2_ROI.label';'L_PreS_ROI.label';'L_ProS_ROI.label';'L_PSL_ROI.label';'L_RI_ROI.label';'L_RSC_ROI.label';'L_s32_ROI.label';'L_s6-8_ROI.label';'L_SCEF_ROI.label';'L_SFL_ROI.label';'L_STGa_ROI.label';'L_STSda_ROI.label';'L_STSdp_ROI.label';'L_STSva_ROI.label';'L_STSvp_ROI.label';'L_STV_ROI.label';'L_TA2_ROI.label';'L_TE1a_ROI.label';'L_TE1m_ROI.label';'L_TE1p_ROI.label';'L_TE2a_ROI.label';'L_TE2p_ROI.label';'L_TF_ROI.label';'L_TGd_ROI.label';'L_TGv_ROI.label';'L_TPOJ1_ROI.label';'L_TPOJ2_ROI.label';'L_TPOJ3_ROI.label';'L_V1_ROI.label';'L_v23ab_ROI.label';'L_V2_ROI.label';'L_V3A_ROI.label';'L_V3B_ROI.label';'L_V3CD_ROI.label';'L_V3_ROI.label';'L_V4_ROI.label';'L_V4t_ROI.label';'L_V6A_ROI.label';'L_V6_ROI.label';'L_V7_ROI.label';'L_V8_ROI.label';'L_VIP_ROI.label';'L_VMV1_ROI.label';'L_VMV2_ROI.label';'L_VMV3_ROI.label';'L_VVC_ROI.label'};
for roi=1:length(llabel_names); 
    try
    %roi=label_names(roi);
    fid = fopen(strcat('/scratch/wbic-beta/pz249/freesurfer_ROCK/preprocessed_fs/', subject_id, '/stats/lh.', llabel_names{roi}, '.stats'));
    linenum = 53;
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    Cstr=C{1, 1}{1, 1};
    fclose(fid);
    D = strsplit(Cstr,' ');
% 1.StructName 2.NumVert 3.SurfArea 4.GrayVol 5.ThickAvg 6.ThickStd 
%%%% 7.MeanCurv 8.GausCurv 9.FoldInd 10.CurvInd

    SurfArea(roi,1)=str2num(D{3});GrayVol(roi,1)=str2num(D{4}); ThickAvg(roi,1)=str2num(D{5});
    MeanCurv(roi,1)=str2num(D{7});GausCurv(roi,1)=str2num(D{8});
    FoldInd(roi,1)=str2num(D{9}); CurvInd(roi,1)=str2num(D{10});   
    end; 
end; %end left loop
Matrix2=horzcat(SurfArea, GrayVol, ThickAvg, MeanCurv, GausCurv, FoldInd, CurvInd);
Matrix2=Matrix2';

% right side
clear SurfArea GrayVol ThickAvg MeanCurv GausCurv FoldInd CurvInd
rlabel_names={'R_10d_ROI.label';'R_10pp_ROI.label';'R_10r_ROI.label';'R_10v_ROI.label';'R_11l_ROI.label';'R_13l_ROI.label';'R_1_ROI.label';'R_23c_ROI.label';'R_23d_ROI.label';'R_24dd_ROI.label';'R_24dv_ROI.label';'R_25_ROI.label';'R_2_ROI.label';'R_31a_ROI.label';'R_31pd_ROI.label';'R_31pv_ROI.label';'R_33pr_ROI.label';'R_3a_ROI.label';'R_3b_ROI.label';'R_43_ROI.label';'R_44_ROI.label';'R_45_ROI.label';'R_46_ROI.label';'R_47l_ROI.label';'R_47m_ROI.label';'R_47s_ROI.label';'R_4_ROI.label';'R_52_ROI.label';'R_55b_ROI.label';'R_5L_ROI.label';'R_5m_ROI.label';'R_5mv_ROI.label';'R_6a_ROI.label';'R_6d_ROI.label';'R_6ma_ROI.label';'R_6mp_ROI.label';'R_6r_ROI.label';'R_6v_ROI.label';'R_7AL_ROI.label';'R_7Am_ROI.label';'R_7m_ROI.label';'R_7PC_ROI.label';'R_7PL_ROI.label';'R_7Pm_ROI.label';'R_8Ad_ROI.label';'R_8Av_ROI.label';'R_8BL_ROI.label';'R_8BM_ROI.label';'R_8C_ROI.label';'R_9-46d_ROI.label';'R_9a_ROI.label';'R_9m_ROI.label';'R_9p_ROI.label';'R_a10p_ROI.label';'R_A1_ROI.label';'R_a24pr_ROI.label';'R_a24_ROI.label';'R_a32pr_ROI.label';'R_a47r_ROI.label';'R_A4_ROI.label';'R_A5_ROI.label';'R_a9-46v_ROI.label';'R_AAIC_ROI.label';'R_AIP_ROI.label';'R_AVI_ROI.label';'R_d23ab_ROI.label';'R_d32_ROI.label';'R_DVT_ROI.label';'R_EC_ROI.label';'R_FEF_ROI.label';'R_FFC_ROI.label';'R_FOP1_ROI.label';'R_FOP2_ROI.label';'R_FOP3_ROI.label';'R_FOP4_ROI.label';'R_FOP5_ROI.label';'R_FST_ROI.label';'R_H_ROI.label';'R_i6-8_ROI.label';'R_IFJa_ROI.label';'R_IFJp_ROI.label';'R_IFSa_ROI.label';'R_IFSp_ROI.label';'R_Ig_ROI.label';'R_IP0_ROI.label';'R_IP1_ROI.label';'R_IP2_ROI.label';'R_IPS1_ROI.label';'R_LBelt_ROI.label';'R_LIPd_ROI.label';'R_LIPv_ROI.label';'R_LO1_ROI.label';'R_LO2_ROI.label';'R_LO3_ROI.label';'R_MBelt_ROI.label';'R_MIP_ROI.label';'R_MI_ROI.label';'R_MST_ROI.label';'R_MT_ROI.label';'R_OFC_ROI.label';'R_OP1_ROI.label';'R_OP2-3_ROI.label';'R_OP4_ROI.label';'R_p10p_ROI.label';'R_p24pr_ROI.label';'R_p24_ROI.label';'R_p32pr_ROI.label';'R_p32_ROI.label';'R_p47r_ROI.label';'R_p9-46v_ROI.label';'R_PBelt_ROI.label';'R_PCV_ROI.label';'R_PeEc_ROI.label';'R_PEF_ROI.label';'R_PFcm_ROI.label';'R_PFm_ROI.label';'R_PFop_ROI.label';'R_PF_ROI.label';'R_PFt_ROI.label';'R_PGi_ROI.label';'R_PGp_ROI.label';'R_PGs_ROI.label';'R_PHA1_ROI.label';'R_PHA2_ROI.label';'R_PHA3_ROI.label';'R_PH_ROI.label';'R_PHT_ROI.label';'R_PI_ROI.label';'R_Pir_ROI.label';'R_PIT_ROI.label';'R_pOFC_ROI.label';'R_PoI1_ROI.label';'R_PoI2_ROI.label';'R_POS1_ROI.label';'R_POS2_ROI.label';'R_PreS_ROI.label';'R_ProS_ROI.label';'R_PSL_ROI.label';'R_RI_ROI.label';'R_RSC_ROI.label';'R_s32_ROI.label';'R_s6-8_ROI.label';'R_SCEF_ROI.label';'R_SFL_ROI.label';'R_STGa_ROI.label';'R_STSda_ROI.label';'R_STSdp_ROI.label';'R_STSva_ROI.label';'R_STSvp_ROI.label';'R_STV_ROI.label';'R_TA2_ROI.label';'R_TE1a_ROI.label';'R_TE1m_ROI.label';'R_TE1p_ROI.label';'R_TE2a_ROI.label';'R_TE2p_ROI.label';'R_TF_ROI.label';'R_TGd_ROI.label';'R_TGv_ROI.label';'R_TPOJ1_ROI.label';'R_TPOJ2_ROI.label';'R_TPOJ3_ROI.label';'R_V1_ROI.label';'R_v23ab_ROI.label';'R_V2_ROI.label';'R_V3A_ROI.label';'R_V3B_ROI.label';'R_V3CD_ROI.label';'R_V3_ROI.label';'R_V4_ROI.label';'R_V4t_ROI.label';'R_V6A_ROI.label';'R_V6_ROI.label';'R_V7_ROI.label';'R_V8_ROI.label';'R_VIP_ROI.label';'R_VMV1_ROI.label';'R_VMV2_ROI.label';'R_VMV3_ROI.label';'R_VVC_ROI.label'};
for roi=1:length(rlabel_names);
    try
    fid = fopen(strcat('/scratch/wbic-beta/pz249/freesurfer_ROCK/preprocessed_fs/', subject_id, '/stats/rh.', rlabel_names{roi}, '.stats'));
    linenum = 53;
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    Cstr=C{1, 1}{1, 1};
    D = strsplit(Cstr,' ');
    fclose(fid);
% 1.StructName 2.NumVert 3.SurfArea 4.GrayVol 5.ThickAvg 6.ThickStd 
%%%% 7.MeanCurv 8.GausCurv 9.FoldInd 10.CurvInd

    SurfArea(roi,1)=str2num(D{3});GrayVol(roi,1)=str2num(D{4}); ThickAvg(roi,1)=str2num(D{5});
    MeanCurv(roi,1)=str2num(D{7});GausCurv(roi,1)=str2num(D{8});
    FoldInd(roi,1)=str2num(D{9}); CurvInd(roi,1)=str2num(D{10});   
    end
end % end right loop
Matrix=horzcat(SurfArea, GrayVol, ThickAvg, MeanCurv, GausCurv, FoldInd, CurvInd);
Matrix=Matrix';

Matrix_full=horzcat(Matrix2, Matrix); % meaning matrix full starts with left and then right!
% z=x-mu_sample/std_sample
for i=1:7
    zscored_matrix(i,:)=(Matrix_full(i,:) - mean(Matrix_full(i,:)))/std(Matrix_full(i,:));
end
R=corr(zscored_matrix);
Rdem=R-mean(mean(R));

Rmat{subject}=Rdem;

MS(subject,:)=mean(Rmat{subject}); %mean similarity

end % end subject loop
label_names_all=vertcat(llabel_names, rlabel_names) %combine all labels

%save('GSK_MS.mat')
%save('Ageing_MS.mat')
%save('IXI_MS_add.mat')
%save('ROCK_MS_add.mat')
%download: rsync -avtry pz249@wbic-gate.vss.cloud.private.cam.ac.uk:/home/pz249/scratch/freesurfer_ATX/preprocessed_fs/ATX_MS.mat /home/brain/Downloads/ATX_MS.mat

%download: rsync -avtry pz249@wbic-gate.vss.cloud.private.cam.ac.uk:/home/pz249/scratch/freesurfer_Ageing/preprocessed_fs/Ageing_MS.mat /home/brain/Downloads/ATX_MS.mat

%download: rsync -avtry pz249@wbic-gate.vss.cloud.private.cam.ac.uk:/home/pz249/scratch/freesurfer_IXI/preprocessed_fs/IXI_MS_add.mat /home/brain/Downloads/IXI_MS_add.mat

%download: rsync -avtry pz249@wbic-gate.vss.cloud.private.cam.ac.uk:/home/pz249/scratch/freesurfer_ROCK/preprocessed_fs/ROCK_MS_add.mat /home/brain/Downloads/IXI_ROCK_add.mat
%% stats on the MS matrices
%dlmwrite('MS.txt', MS(:, 181:360))
mean_MS=mean(MS) % group mean similarity of each ROI!
%ATX
for roi=1:360
    [h,p, ci, stats]=ttest2(MS(1:39, roi), MS(40:79, roi));
    tstat_roi(roi)=stats.tstat;
    pstat_roi(roi)=p;
end

%Ageing
for roi=1:360
    [h,p, ci, stats]=ttest2(MS(1:31, roi), MS(32:65, roi));
    tstat_roi(roi)=stats.tstat;
    pstat_roi(roi)=p;
end



FDR = mafdr(pstat_roi(181:360))
sort(FDR)
sort(pstat_roi(181:360))

min(tstat_roi); max(tstat_roi)
label_names_all(tstat_roi>2.5)
label_names_all(tstat_roi<-2.2)
label_names_all(pstat_roi<0.01)
mean(MS(1:49, 262))
mean(MS(50:95, 262))
max(mean_MS)
label_names_all(mean_MS>0.04)
label_names_all(mean_MS<-0.06)


%plot(mean(Rdem))
%imagesc(Rdem); save('MS.mat')
%figure; imagesc(Rmat{1})

%mean MS across all ROIs:
totalMS=mean(MS');
[h,p]=ttest2(totalMS(1:39), totalMS(40:79)) % => Not Significant


%% first results:

    %'L_TPOJ1_ROI.label'  - temporo-parietal occipital junction
    %'R_6v_ROI.label' - precentral gyrus ventral
    %'R_AVI_ROI.label' - Insula
    %'R_MI_ROI.label' - Insula
    %'R_p32pr_ROI.label' - ACC
    
    %'R_10r_ROI.label' - aPFC
    %'R_IFSa_ROI.label' - IFG
    %'R_IFSp_ROI.label' - IFG

    
%% COMBINED
label_names_all=vertcat(llabel_names, rlabel_names) %combine all labels
%save('All_MS.mat')

mean_MS=mean(MS) % group mean similarity of each ROI!
%ATX
for roi=1:360
    [h,p, ci, stats]=ttest2(MS(1:148, roi), MS(149:331, roi));
    tstat_roi(roi)=stats.tstat;
    pstat_roi(roi)=p;
end

min(tstat_roi); max(tstat_roi)
label_names_all(tstat_roi>2.5)
label_names_all(tstat_roi<-2.5)
label_names_all(pstat_roi<0.01)
mean(MS(1:49, 262))
mean(MS(50:95, 262))
max(mean_MS)
%!!!
FDR = mafdr(pstat_roi(1:360))
label_names_all(FDR<0.01)
rlabel_names(FDR<0.05)

%
   % 'L_25_ROI.label'
    %'L_3a_ROI.label'
    %'L_3b_ROI.label'
    %'L_47s_ROI.label'
    %'L_a10p_ROI.label'
    %'L_AIP_ROI.label'
    %'L_AVI_ROI.label'
    %'L_d23ab_ROI.label'
    %'L_IP1_ROI.label'
    %'L_OFC_ROI.label'
    %'L_POS1_ROI.label'
    %'L_POS2_ROI.label'
    %'L_STGa_ROI.label'
    %'L_V3A_ROI.label'
    %'L_V6_ROI.label'
    %'R_10r_ROI.label'
    %'R_3a_ROI.label'
    %'R_3b_ROI.label'
    %'R_45_ROI.label'
    %'R_6mp_ROI.label'
    %'R_6v_ROI.label'
    %'R_7AL_ROI.label'
    %'R_9-46d_ROI.label'
    %'R_a10p_ROI.label'
    %'R_A5_ROI.label'
    %'R_a9-46v_ROI.label'
    %'R_AAIC_ROI.label'
    %'R_AIP_ROI.label'
    %'R_AVI_ROI.label'
    %'R_IFSp_ROI.label'
    %%'R_MI_ROI.label'
    %'R_OFC_ROI.label'
    %'R_p10p_ROI.label'
    %'R_p32_ROI.label'
    %'R_p9-46v_ROI.label'
    %'R_PeEc_ROI.label'
    %'R_PHT_ROI.label'
    %'R_POS2_ROI.label'
    %'R_SFL_ROI.label'
    %'R_STGa_ROI.label'
    %'R_TE1m_ROI.label'
    %'R_v23ab_ROI.label'
    %'R_V4t_ROI.label'%
    
    
