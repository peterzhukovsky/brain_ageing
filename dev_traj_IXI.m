%% %% 1. constructing developmental trajectories for each ROI across all ages
%No. Label Name:				R   G   B   A
llabel_names={'L_10d_ROI.label';'L_10pp_ROI.label';'L_10r_ROI.label';'L_10v_ROI.label';'L_11l_ROI.label';'L_13l_ROI.label';'L_1_ROI.label';'L_23c_ROI.label';'L_23d_ROI.label';'L_24dd_ROI.label';'L_24dv_ROI.label';'L_25_ROI.label';'L_2_ROI.label';'L_31a_ROI.label';'L_31pd_ROI.label';'L_31pv_ROI.label';'L_33pr_ROI.label';'L_3a_ROI.label';'L_3b_ROI.label';'L_43_ROI.label';'L_44_ROI.label';'L_45_ROI.label';'L_46_ROI.label';'L_47l_ROI.label';'L_47m_ROI.label';'L_47s_ROI.label';'L_4_ROI.label';'L_52_ROI.label';'L_55b_ROI.label';'L_5L_ROI.label';'L_5m_ROI.label';'L_5mv_ROI.label';'L_6a_ROI.label';'L_6d_ROI.label';'L_6ma_ROI.label';'L_6mp_ROI.label';'L_6r_ROI.label';'L_6v_ROI.label';'L_7AL_ROI.label';'L_7Am_ROI.label';'L_7m_ROI.label';'L_7PC_ROI.label';'L_7PL_ROI.label';'L_7Pm_ROI.label';'L_8Ad_ROI.label';'L_8Av_ROI.label';'L_8BL_ROI.label';'L_8BM_ROI.label';'L_8C_ROI.label';'L_9-46d_ROI.label';'L_9a_ROI.label';'L_9m_ROI.label';'L_9p_ROI.label';'L_a10p_ROI.label';'L_A1_ROI.label';'L_a24pr_ROI.label';'L_a24_ROI.label';'L_a32pr_ROI.label';'L_a47r_ROI.label';'L_A4_ROI.label';'L_A5_ROI.label';'L_a9-46v_ROI.label';'L_AAIC_ROI.label';'L_AIP_ROI.label';'L_AVI_ROI.label';'L_d23ab_ROI.label';'L_d32_ROI.label';'L_DVT_ROI.label';'L_EC_ROI.label';'L_FEF_ROI.label';'L_FFC_ROI.label';'L_FOP1_ROI.label';'L_FOP2_ROI.label';'L_FOP3_ROI.label';'L_FOP4_ROI.label';'L_FOP5_ROI.label';'L_FST_ROI.label';'L_H_ROI.label';'L_i6-8_ROI.label';'L_IFJa_ROI.label';'L_IFJp_ROI.label';'L_IFSa_ROI.label';'L_IFSp_ROI.label';'L_Ig_ROI.label';'L_IP0_ROI.label';'L_IP1_ROI.label';'L_IP2_ROI.label';'L_IPS1_ROI.label';'L_LBelt_ROI.label';'L_LIPd_ROI.label';'L_LIPv_ROI.label';'L_LO1_ROI.label';'L_LO2_ROI.label';'L_LO3_ROI.label';'L_MBelt_ROI.label';'L_MIP_ROI.label';'L_MI_ROI.label';'L_MST_ROI.label';'L_MT_ROI.label';'L_OFC_ROI.label';'L_OP1_ROI.label';'L_OP2-3_ROI.label';'L_OP4_ROI.label';'L_p10p_ROI.label';'L_p24pr_ROI.label';'L_p24_ROI.label';'L_p32pr_ROI.label';'L_p32_ROI.label';'L_p47r_ROI.label';'L_p9-46v_ROI.label';'L_PBelt_ROI.label';'L_PCV_ROI.label';'L_PeEc_ROI.label';'L_PEF_ROI.label';'L_PFcm_ROI.label';'L_PFm_ROI.label';'L_PFop_ROI.label';'L_PF_ROI.label';'L_PFt_ROI.label';'L_PGi_ROI.label';'L_PGp_ROI.label';'L_PGs_ROI.label';'L_PHA1_ROI.label';'L_PHA2_ROI.label';'L_PHA3_ROI.label';'L_PH_ROI.label';'L_PHT_ROI.label';'L_PI_ROI.label';'L_Pir_ROI.label';'L_PIT_ROI.label';'L_pOFC_ROI.label';'L_PoI1_ROI.label';'L_PoI2_ROI.label';'L_POS1_ROI.label';'L_POS2_ROI.label';'L_PreS_ROI.label';'L_ProS_ROI.label';'L_PSL_ROI.label';'L_RI_ROI.label';'L_RSC_ROI.label';'L_s32_ROI.label';'L_s6-8_ROI.label';'L_SCEF_ROI.label';'L_SFL_ROI.label';'L_STGa_ROI.label';'L_STSda_ROI.label';'L_STSdp_ROI.label';'L_STSva_ROI.label';'L_STSvp_ROI.label';'L_STV_ROI.label';'L_TA2_ROI.label';'L_TE1a_ROI.label';'L_TE1m_ROI.label';'L_TE1p_ROI.label';'L_TE2a_ROI.label';'L_TE2p_ROI.label';'L_TF_ROI.label';'L_TGd_ROI.label';'L_TGv_ROI.label';'L_TPOJ1_ROI.label';'L_TPOJ2_ROI.label';'L_TPOJ3_ROI.label';'L_V1_ROI.label';'L_v23ab_ROI.label';'L_V2_ROI.label';'L_V3A_ROI.label';'L_V3B_ROI.label';'L_V3CD_ROI.label';'L_V3_ROI.label';'L_V4_ROI.label';'L_V4t_ROI.label';'L_V6A_ROI.label';'L_V6_ROI.label';'L_V7_ROI.label';'L_V8_ROI.label';'L_VIP_ROI.label';'L_VMV1_ROI.label';'L_VMV2_ROI.label';'L_VMV3_ROI.label';'L_VVC_ROI.label'};
rlabel_names={'R_10d_ROI.label';'R_10pp_ROI.label';'R_10r_ROI.label';'R_10v_ROI.label';'R_11l_ROI.label';'R_13l_ROI.label';'R_1_ROI.label';'R_23c_ROI.label';'R_23d_ROI.label';'R_24dd_ROI.label';'R_24dv_ROI.label';'R_25_ROI.label';'R_2_ROI.label';'R_31a_ROI.label';'R_31pd_ROI.label';'R_31pv_ROI.label';'R_33pr_ROI.label';'R_3a_ROI.label';'R_3b_ROI.label';'R_43_ROI.label';'R_44_ROI.label';'R_45_ROI.label';'R_46_ROI.label';'R_47l_ROI.label';'R_47m_ROI.label';'R_47s_ROI.label';'R_4_ROI.label';'R_52_ROI.label';'R_55b_ROI.label';'R_5L_ROI.label';'R_5m_ROI.label';'R_5mv_ROI.label';'R_6a_ROI.label';'R_6d_ROI.label';'R_6ma_ROI.label';'R_6mp_ROI.label';'R_6r_ROI.label';'R_6v_ROI.label';'R_7AL_ROI.label';'R_7Am_ROI.label';'R_7m_ROI.label';'R_7PC_ROI.label';'R_7PL_ROI.label';'R_7Pm_ROI.label';'R_8Ad_ROI.label';'R_8Av_ROI.label';'R_8BL_ROI.label';'R_8BM_ROI.label';'R_8C_ROI.label';'R_9-46d_ROI.label';'R_9a_ROI.label';'R_9m_ROI.label';'R_9p_ROI.label';'R_a10p_ROI.label';'R_A1_ROI.label';'R_a24pr_ROI.label';'R_a24_ROI.label';'R_a32pr_ROI.label';'R_a47r_ROI.label';'R_A4_ROI.label';'R_A5_ROI.label';'R_a9-46v_ROI.label';'R_AAIC_ROI.label';'R_AIP_ROI.label';'R_AVI_ROI.label';'R_d23ab_ROI.label';'R_d32_ROI.label';'R_DVT_ROI.label';'R_EC_ROI.label';'R_FEF_ROI.label';'R_FFC_ROI.label';'R_FOP1_ROI.label';'R_FOP2_ROI.label';'R_FOP3_ROI.label';'R_FOP4_ROI.label';'R_FOP5_ROI.label';'R_FST_ROI.label';'R_H_ROI.label';'R_i6-8_ROI.label';'R_IFJa_ROI.label';'R_IFJp_ROI.label';'R_IFSa_ROI.label';'R_IFSp_ROI.label';'R_Ig_ROI.label';'R_IP0_ROI.label';'R_IP1_ROI.label';'R_IP2_ROI.label';'R_IPS1_ROI.label';'R_LBelt_ROI.label';'R_LIPd_ROI.label';'R_LIPv_ROI.label';'R_LO1_ROI.label';'R_LO2_ROI.label';'R_LO3_ROI.label';'R_MBelt_ROI.label';'R_MIP_ROI.label';'R_MI_ROI.label';'R_MST_ROI.label';'R_MT_ROI.label';'R_OFC_ROI.label';'R_OP1_ROI.label';'R_OP2-3_ROI.label';'R_OP4_ROI.label';'R_p10p_ROI.label';'R_p24pr_ROI.label';'R_p24_ROI.label';'R_p32pr_ROI.label';'R_p32_ROI.label';'R_p47r_ROI.label';'R_p9-46v_ROI.label';'R_PBelt_ROI.label';'R_PCV_ROI.label';'R_PeEc_ROI.label';'R_PEF_ROI.label';'R_PFcm_ROI.label';'R_PFm_ROI.label';'R_PFop_ROI.label';'R_PF_ROI.label';'R_PFt_ROI.label';'R_PGi_ROI.label';'R_PGp_ROI.label';'R_PGs_ROI.label';'R_PHA1_ROI.label';'R_PHA2_ROI.label';'R_PHA3_ROI.label';'R_PH_ROI.label';'R_PHT_ROI.label';'R_PI_ROI.label';'R_Pir_ROI.label';'R_PIT_ROI.label';'R_pOFC_ROI.label';'R_PoI1_ROI.label';'R_PoI2_ROI.label';'R_POS1_ROI.label';'R_POS2_ROI.label';'R_PreS_ROI.label';'R_ProS_ROI.label';'R_PSL_ROI.label';'R_RI_ROI.label';'R_RSC_ROI.label';'R_s32_ROI.label';'R_s6-8_ROI.label';'R_SCEF_ROI.label';'R_SFL_ROI.label';'R_STGa_ROI.label';'R_STSda_ROI.label';'R_STSdp_ROI.label';'R_STSva_ROI.label';'R_STSvp_ROI.label';'R_STV_ROI.label';'R_TA2_ROI.label';'R_TE1a_ROI.label';'R_TE1m_ROI.label';'R_TE1p_ROI.label';'R_TE2a_ROI.label';'R_TE2p_ROI.label';'R_TF_ROI.label';'R_TGd_ROI.label';'R_TGv_ROI.label';'R_TPOJ1_ROI.label';'R_TPOJ2_ROI.label';'R_TPOJ3_ROI.label';'R_V1_ROI.label';'R_v23ab_ROI.label';'R_V2_ROI.label';'R_V3A_ROI.label';'R_V3B_ROI.label';'R_V3CD_ROI.label';'R_V3_ROI.label';'R_V4_ROI.label';'R_V4t_ROI.label';'R_V6A_ROI.label';'R_V6_ROI.label';'R_V7_ROI.label';'R_V8_ROI.label';'R_VIP_ROI.label';'R_VMV1_ROI.label';'R_VMV2_ROI.label';'R_VMV3_ROI.label';'R_VVC_ROI.label'};
label_names_all=vertcat(llabel_names, rlabel_names) %combine all labels
label_names_figure={'left 10d';'left 10pp';'left 10';'left 10v';'left 11';'left 13';'left 1';'left 23c';'left 23d';'left 24dd';'left 24dv';'left 25';'left 2';'left 31a';'left 31pd';'left 31pv';'left 33p';'left 3a';'left 3b';'left 43';'left 44';'left 45';'left 46';'left 47';'left 47m';'left 47s';'left 4';'left 52';'left 55b';'left 5';'left 5m';'left 5mv';'left 6a';'left 6d';'left 6ma';'left 6mp';'left 6';'left 6v';'left 7A';'left 7Am';'left 7m';'left 7PC';'left 7P';'left 7Pm';'left 8Ad';'left 8Av';'left 8B';'left 8BM';'left 8C';'left 9-46d';'left 9a';'left 9m';'left 9p';'left a10p';'left A1';'left a24p';'left a24';'left a32p';'left a47';'left A4';'left A5';'left a9-46v';'left AAIC';'left AIP';'left AVI';'left d23ab';'left d32';'left DVT';'left EC';'left FEF';'left FFC';'left FOP1';'left FOP2';'left FOP3';'left FOP4';'left FOP5';'left FST';'left H';'left i6-8';'left IFJa';'left IFJp';'left IFSa';'left IFSp';'left Ig';'left IP0';'left IP1';'left IP2';'left IPS1';'left LBelt';'left LIPd';'left LIPv';'left LO1';'left LO2';'left LO3';'left MBelt';'left MIP';'left MI';'left MST';'left MT';'left OFC';'left OP1';'left OP2-3';'left OP4';'left p10p';'left p24p';'left p24';'left p32p';'left p32';'left p47';'left p9-46v';'left PBelt';'left PCV';'left PeEc';'left PEF';'left PFcm';'left PFm';'left PFop';'left PF';'left PFt';'left PGi';'left PGp';'left PGs';'left PHA1';'left PHA2';'left PHA3';'left PH';'left PHT';'left PI';'left Pi';'left PIT';'left pOFC';'left PoI1';'left PoI2';'left POS1';'left POS2';'left PreS';'left ProS';'left PS';'left RI';'left RSC';'left s32';'left s6-8';'left SCEF';'left SF';'left STGa';'left STSda';'left STSdp';'left STSva';'left STSvp';'left STV';'left TA2';'left TE1a';'left TE1m';'left TE1p';'left TE2a';'left TE2p';'left TF';'left TGd';'left TGv';'left TPOJ1';'left TPOJ2';'left TPOJ3';'left V1';'left v23ab';'left V2';'left V3A';'left V3B';'left V3CD';'left V3';'left V4';'left V4t';'left V6A';'left V6';'left V7';'left V8';'left VIP';'left VMV1';'left VMV2';'left VMV3';'left VVC';'right 10d';'right 10pp';'right 10';'right 10v';'right 11';'right 13';'right 1';'right 23c';'right 23d';'right 24dd';'right 24dv';'right 25';'right 2';'right 31a';'right 31pd';'right 31pv';'right 33p';'right 3a';'right 3b';'right 43';'right 44';'right 45';'right 46';'right 47';'right 47m';'right 47s';'right 4';'right 52';'right 55b';'right 5';'right 5m';'right 5mv';'right 6a';'right 6d';'right 6ma';'right 6mp';'right 6';'right 6v';'right 7A';'right 7Am';'right 7m';'right 7PC';'right 7P';'right 7Pm';'right 8Ad';'right 8Av';'right 8B';'right 8BM';'right 8C';'right 9-46d';'right 9a';'right 9m';'right 9p';'right a10p';'right A1';'right a24p';'right a24';'right a32p';'right a47';'right A4';'right A5';'right a9-46v';'right AAIC';'right AIP';'right AVI';'right d23ab';'right d32';'right DVT';'right EC';'right FEF';'right FFC';'right FOP1';'right FOP2';'right FOP3';'right FOP4';'right FOP5';'right FST';'right H';'right i6-8';'right IFJa';'right IFJp';'right IFSa';'right IFSp';'right Ig';'right IP0';'right IP1';'right IP2';'right IPS1';'right LBelt';'right LIPd';'right LIPv';'right LO1';'right LO2';'right LO3';'right MBelt';'right MIP';'right MI';'right MST';'right MT';'right OFC';'right OP1';'right OP2-3';'right OP4';'right p10p';'right p24p';'right p24';'right p32p';'right p32';'right p47';'right p9-46v';'right PBelt';'right PCV';'right PeEc';'right PEF';'right PFcm';'right PFm';'right PFop';'right PF';'right PFt';'right PGi';'right PGp';'right PGs';'right PHA1';'right PHA2';'right PHA3';'right PH';'right PHT';'right PI';'right Pi';'right PIT';'right pOFC';'right PoI1';'right PoI2';'right POS1';'right POS2';'right PreS';'right ProS';'right PS';'right RI';'right RSC';'right s32';'right s6-8';'right SCEF';'right SF';'right STGa';'right STSda';'right STSdp';'right STSva';'right STSvp';'right STV';'right TA2';'right TE1a';'right TE1m';'right TE1p';'right TE2a';'right TE2p';'right TF';'right TGd';'right TGv';'right TPOJ1';'right TPOJ2';'right TPOJ3';'right V1';'right v23ab';'right V2';'right V3A';'right V3B';'right V3CD';'right V3';'right V4';'right V4t';'right V6A';'right V6';'right V7';'right V8';'right VIP';'right VMV1';'right VMV2';'right VMV3';'right VVC'};
MS=zeros(1, 148); age=zeros(1,1); group=zeros(1,1); datasetNr=0; yrs_use=0; sex=0;
MS(~(age>0),:)=[];group(~(age>0))=[]; age(~(age>0))=[];id=1:length(age);
%https://uk.mathworks.com/help/matlab/function-handles.html

for roi=1:360%
    %subplot(5,6,roi)
    %scatter(age,MS(1:148,roi))
    %https://uk.mathworks.com/help/curvefit/fit.html [fitobject,gof] = fit(x,y,fitType)
    x=age(group==4); y=MS(group==4,roi);
    [p,S]=polyfit(x,y,2); p_all(:,roi)=p;
    %f = @(x) p(1)*x^2+p(2)*x; % + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    %fplot(f);hold on; %scatter(age,MS(:,roi));  %xlim([20, 70]); ylim([-0.05 0.035]); 
    
    [poly1,gof1] = fit(x,y,'poly1');
    [poly2,gof2] = fit(x,y,'poly2');
    [poly3,gof3] = fit(x,y,'poly3');
    adjr1(roi)=gof1.adjrsquare; r1(roi)=gof1.rsquare; 
    adjr2(roi)=gof2.adjrsquare; r2(roi)=gof2.rsquare; 
    adjr3(roi)=gof3.adjrsquare; r3(roi)=gof3.rsquare; 
end


%% hierarchical clustering - not v useful; grouping the trajectories by U vs flat vs inverted U
mean(adjr2(FDR<0.01))
mean(r2(FDR<0.01))
figure; histogram(adjr2)    
label_names_all(adjr2<0 & FDR<0.01)
label_names_all(FDR<0.05)

scatter((p_all(1,:)), p_all(2,:))
C_data=vertcat(log(p_all(1,:)), p_all(2,:))

Y = pdist(cluster_val);%Y=pdist(p_all(1:2,:)');% 
squareform(Y); Z = linkage(Y); figure; dendrogram(Z)
T = cluster(Z,'cutoff',1.155); T = cluster(Z,'maxclust',6)

%%%%% manual clustering by differentiating at 4 age points 
for roi=1:360
p=p_all(:,roi);
f = @(x) p(1)*x^2+p(2)*x + p(3);
df = @(x) p(1)*2*x+p(2); x0=-p(2)/p(1)/2; dfx0(roi)=x0;
cluster_val(roi, 1)=feval(df, 27); 
cluster_val(roi, 2)=feval(df, 35);
cluster_val(roi, 3)=feval(df, 45);
cluster_val(roi, 4)=feval(df, 60);
end
cluster_val=cluster_val*10^4;
%cluster_val=cluster_val(FDR<0.01,:)
figure;scatter(cluster_val(:,1),cluster_val(:,4))
T(cluster_val(:,1)>0 & cluster_val(:,4)>0)=1
T(cluster_val(:,1)<0 & cluster_val(:,4)<0)=2
T(cluster_val(:,1)>0 & cluster_val(:,4)<0)=3
T(cluster_val(:,1)<0 & cluster_val(:,4)>0)=4

%https://uk.mathworks.com/help/stats/hierarchical-clustering.html

%% plotting the trajectories
id=1:360';
figure(2); hold off; %swap between plotting all lines in grey and plotting colored lines over them
for r=1:length(T);
    roi=id(FDR<0.05); %243 1:360;% 
    roi=roi(r);
    p=p_all(:,roi);
    f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    %if T(roi)==1, COL = {'b'}; elseif T(roi)==2, COL = {'-r'}; elseif T(roi)==3, COL = {'-k'} ; elseif T(roi)==4, COL = {'-g'}; end;
    COL={'b'};
    fp=fplot(f, [20 80], COL{:});hold on; % h=fplot(f, [20 80]); set(h, {'color'}, {[T(r)/20 T(r)/10 T(r)/25]}); hold on; %
    fp.Color = [0.85 0.85 0.85];
end
for r=1:length(T);
    roi=id(T==2 & FDR<0.05); %243 1:360;% 
    roi=roi(r);
    p=p_all(:,roi);
    f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    if T(roi)==1, COL = {'b'}; elseif T(roi)==2, COL = {'-r'}; elseif T(roi)==3, COL = {'-k'} ; elseif T(roi)==4, COL = {'-g'}; end;
    fp=fplot(f, [20 80], COL{:});hold on; % h=fplot(f, [20 80]); set(h, {'color'}, {[T(r)/20 T(r)/10 T(r)/25]}); hold on; %
    %fp.Color = [0.85 0.85 0.85];
end
 
%% Plotting one ROI for method viz
figure(3); hold off;
for roi=243
   p=p_all(:,roi);
    f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    %if T(r)==1, COL = {'-b'}; elseif T(r)==2, COL = {'-r'}; elseif T(r)==3, COL = {'-k'} ; elseif T(r)==4, COL = {'-g'}, end;
    
    fplot(f, [20 80]);hold on; % h=fplot(f, [20 80]); set(h, {'color'}, {[T(r)/20 T(r)/10 T(r)/25]}); hold on; %
        scatter(age(group==3),MS(group==3,roi), 20, 'k');  
        x=age(group==3); CI = predint(poly2,x,0.95,'functional','on'); plot(x,CI,'k.', 'MarkerSize',3); xlim([18 80]);
        scatter(age(group==2),MS((group==2),roi), 10, 'r', 'filled'); hold on;
        scatter(age(group==1),MS((group==1),roi), 10, 'b', 'filled'); 
    
end

%% MS of each participant compared to normative development trajectories and put into categories
id=1:964;
for subject=(id(group<=2));
    for roi=1:360
    p=p_all(:,roi);
    f = @(x) p(1)*x^2+p(2)*x + p(3);
    MS_pred(subject, roi) = feval(f,age(subject));
    end
end

MS_difference=MS(group==2,:)- MS_pred(group==2, :);
MS_difference_HC=MS(group==1,:)- MS_pred(group==1, :);

mean_MS_difference=mean(MS_difference);
mean_MS_difference_HC=mean(MS_difference_HC); figure(12); histogram(mean_MS_difference, 25)
mean_MS_difference(FDR<0.01)

age_SDI=age(group==2); age_HC=age(group==1);

%put the MS_diff variable into categories - 8 categories are too fine for the 180 participants we have
for r=1:360% length(T)% 
    roi=1:360;% id(FDR<0.01); %243 %
    roi=roi(r);
MS_diff_cat_age(r,1)=mean(MS_difference(age_SDI>18 & age_SDI<=30,roi)); %mean age=25.514 mean(age_SDI(age_SDI>18 & age_SDI<=30))
MS_diff_cat_age(r,2)=mean(MS_difference(age_SDI>30 & age_SDI<=35,roi)); %mean=33.0638 mean(age_SDI(age_SDI>30 & age_SDI<=35))
MS_diff_cat_age(r,3)=mean(MS_difference(age_SDI>35 & age_SDI<=40,roi)); %mean=37.7058 mean(age_SDI(age_SDI>35 & age_SDI<=40))
MS_diff_cat_age(r,4)=mean(MS_difference(age_SDI>40 & age_SDI<=45,roi)); %mean=42.8025 mean(age_SDI(age_SDI>40 & age_SDI<=45))
MS_diff_cat_age(r,5)=mean(MS_difference(age_SDI>45 & age_SDI<=60,roi)); %mean=49.1819 mean(age_SDI(age_SDI>45 & age_SDI<=60))

MS_diff_cat_HC(r,1)=mean(MS_difference_HC(age_HC>18 & age_HC<=30,roi)); %mean age=25.2145 mean(age_HC(age_HC>18 & age_HC<=30))
MS_diff_cat_HC(r,2)=mean(MS_difference_HC(age_HC>30 & age_HC<=35,roi)); %mean=35.6275 mean(age_HC(age_HC>30 & age_HC<=40))
MS_diff_cat_HC(r,3)=mean(MS_difference_HC(age_HC>35 & age_HC<=40,roi)); %mean=47.3538 mean(age_HC(age_HC>40 & age_HC<=60))
MS_diff_cat_HC(r,4)=mean(MS_difference_HC(age_HC>40 & age_HC<=45,roi));
MS_diff_cat_HC(r,5)=mean(MS_difference_HC(age_HC>45 & age_HC<=60,roi));
end

histogram(T(FDR<0.01))

figure(2); hold off
for r=1:20%length(T)
    roi=id(FDR<0.01 & T==4); %1:360;%id(FDR<0.01); %243 %
    roi=roi(r);
    subplot(4,1,r)
    %if FDR(roi)<0.01, COL = 'r'; elseif FDR(roi)>=0.01, COL = 'b'; end; %elseif T(r)==3, COL = {'-k'} ; elseif T(r)==4, COL = {'-g'},
    COL = 'r'; bar(MS_diff_cat_age(roi,:), COL); set(gca,'xtick',[]); xlabel(label_names_figure(roi)); hold on; ylim([-0.05 0.05]); set(gca,'YTickLabel',[]); 
    COL = 'b'; b=bar(MS_diff_cat_HC(roi,:), COL); set(b,'FaceAlpha',.55); set(gca,'xtick',[]); xlabel(label_names_figure(roi)); hold on; ylim([-0.05 0.05]); set(gca,'YTickLabel',[]); 
    p=p_all(:,roi);
    f = @(x) p(1)*x^2+p(2)*x ;%+ p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
    ax1 = gca; % current axes 
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
    hold on; fplot(f, [17 70], '-b', 'LineWidth', 3); set(gca,'YTickLabel',[]); xlim([18 68]); set(gca,'xtick',[]); hold off; %
end






%% t-testing the MS_difference
% SDI
for roi=1:360
    [h,p,ci,stats] = ttest(MS_difference(age_SDI>18 & age_SDI<=30,roi));
    tstat1(roi)=stats.tstat; pstat1(roi)=p; %%%1
    [h,p,ci,stats] = ttest(MS_difference(age_SDI>30 & age_SDI<=35,roi));
    tstat2(roi)=stats.tstat; pstat2(roi)=p;%%%2
    [h,p,ci,stats] = ttest(MS_difference(age_SDI>35 & age_SDI<=40,roi));
    tstat3(roi)=stats.tstat; pstat3(roi)=p;%%%3
    [h,p,ci,stats] = ttest(MS_difference(age_SDI>40 & age_SDI<=45,roi));
    tstat4(roi)=stats.tstat; pstat4(roi)=p;%%%4
    [h,p,ci,stats] = ttest(MS_difference(age_SDI>45 & age_SDI<=60,roi));
    tstat5(roi)=stats.tstat; pstat5(roi)=p;%%%5
end
FDR1=mafdr(pstat1); FDR2=mafdr(pstat2); FDR3=mafdr(pstat3);FDR4=mafdr(pstat4); FDR5=mafdr(pstat5);
    %%% consistent tstats are those that a) are FDR corrected and b) in the
    %%% same direction as the T class
consistent1=label_names_all((tstat1>0 & (T==1 |T==3) & FDR1<0.01) | (tstat1<0 & (T==2 |T==4) & FDR1<0.01)& adjr2>0) %c 
inconsistent1=label_names_all((tstat1<0 & (T==1 |T==3) & FDR1<0.01) | (tstat1>0 & (T==2 |T==4) & FDR1<0.01)& adjr2>0) %inc
consistent2=label_names_all((tstat2>0 & (T==1 |T==4) & FDR2<0.01) | (tstat2<0 & (T==2 |T==3) & FDR2<0.01)& adjr2>0)
inconsistent2=label_names_all((tstat2<0 & (T==1 |T==4) & FDR2<0.01) | (tstat2>0 & (T==2 |T==3) & FDR2<0.01)& adjr2>0) %inc
consistent3=label_names_all((tstat3>0 & (T==1 |T==4) & FDR3<0.01) | (tstat3<0 & (T==2 |T==3) & FDR3<0.01)& adjr2>0)
inconsistent3=label_names_all((tstat3<0 & (T==1 |T==4) & FDR3<0.01) | (tstat3>0 & (T==2 |T==3) & FDR3<0.01)& adjr2>0) %inc
consistent4=label_names_all((tstat4>0 & (T==1 |T==4) & FDR4<0.01) | (tstat4<0 & (T==2 |T==3) & FDR4<0.01)& adjr2>0)
inconsistent4=label_names_all((tstat4<0 & (T==1 |T==4) & FDR4<0.01) | (tstat4>0 & (T==2 |T==3) & FDR4<0.01)& adjr2>0) %inc
consistent5=label_names_all((tstat5>0 & (T==1 |T==4) & FDR5<0.01) | (tstat5<0 & (T==2 |T==3) & FDR5<0.01)& adjr2>0)
inconsistent5=label_names_all((tstat5<0 & (T==1 |T==4) & FDR5<0.01) | (tstat5>0 & (T==2 |T==3) & FDR5<0.01)& adjr2>0) %inc

consistent1=((tstat1>0 & (T==1 |T==3) & FDR1<0.01) | (tstat1<0 & (T==2 |T==4) & FDR1<0.01)) %c 
inconsistent1=((tstat1<0 & (T==1 |T==3) & FDR1<0.01) | (tstat1>0 & (T==2 |T==4) & FDR1<0.01)) %inc
consistent2=((tstat2>0 & (T==1 |T==4) & FDR2<0.01) | (tstat2<0 & (T==2 |T==3) & FDR2<0.01))
inconsistent2=((tstat2<0 & (T==1 |T==4) & FDR2<0.01) | (tstat2>0 & (T==2 |T==3) & FDR2<0.01)) %inc
consistent3=((tstat3>0 & (T==1 |T==4) & FDR3<0.01) | (tstat3<0 & (T==2 |T==3) & FDR3<0.01))
inconsistent3=((tstat3<0 & (T==1 |T==4) & FDR3<0.01) | (tstat3>0 & (T==2 |T==3) & FDR3<0.01)) %inc
consistent4=((tstat4>0 & (T==1 |T==4) & FDR4<0.01) | (tstat4<0 & (T==2 |T==3) & FDR4<0.01))
inconsistent4=((tstat4<0 & (T==1 |T==4) & FDR4<0.01) | (tstat4>0 & (T==2 |T==3) & FDR4<0.01)) %inc
consistent5=((tstat5>0 & (T==1 |T==4) & FDR5<0.01) | (tstat5<0 & (T==2 |T==3) & FDR5<0.01))
inconsistent5=((tstat5<0 & (T==1 |T==4) & FDR5<0.01) | (tstat5>0 & (T==2 |T==3) & FDR5<0.01)) %inc


% HC
    for roi=1:360
        [h,p,ci,stats] = ttest(MS_difference_HC(age_HC>18 & age_HC<=30,roi));
        tstat1(roi)=stats.tstat; pstat1(roi)=p; %%%1
        [h,p,ci,stats] = ttest(MS_difference_HC(age_HC>30 & age_HC<=35,roi));
        tstat2(roi)=stats.tstat; pstat2(roi)=p;%%%2
        [h,p,ci,stats] = ttest(MS_difference_HC(age_HC>35 & age_HC<=40,roi));
        tstat3(roi)=stats.tstat; pstat3(roi)=p;%%%3
        [h,p,ci,stats] = ttest(MS_difference_HC(age_HC>40 & age_HC<=45,roi));
        tstat4(roi)=stats.tstat; pstat4(roi)=p;%%%4
        [h,p,ci,stats] = ttest(MS_difference_HC(age_HC>45 & age_HC<=60,roi));
        tstat5(roi)=stats.tstat; pstat5(roi)=p;%%%5
    end

FDR1=mafdr(pstat1); FDR2=mafdr(pstat2); FDR3=mafdr(pstat3); FDR4=mafdr(pstat4); FDR5=mafdr(pstat5);





%% Prediction analysis (not enough data and relationships are nonlinear)

%for each 5 year interval see which age the function takes you to based on
%the MS_diff_roi

for r=1:length(id(FDR<0.01))
    roi=id(FDR<0.01); %243%1:360
    roi=roi(r);
    
    p=p_all(:,roi); f = @(x) p(1)*x^2+p(2)*x + p(3);  
    p_new= feval(f,25.514) + MS_diff_cat_age(roi, 1);
    p(3)=p(3)-p_new;
  	try
        if max(roots(p))<150, ageing_metric(r, 1)=max(roots(p)); elseif min(roots(p))>30 ageing_metric(r, 1)=min(roots(p)); end; 
    end
    
    p=p_all(:,roi); f = @(x) p(1)*x^2+p(2)*x + p(3);  
    p_new= feval(f,35.4157) + MS_diff_cat_age(roi, 1);
    p(3)=p(3)-p_new;
    try
        if max(roots(p))<150, ageing_metric(r, 2)=max(roots(p)); elseif min(roots(p))>30 ageing_metric(r, 2)=min(roots(p)); end; 
    end
    
    p=p_all(:,roi); f = @(x) p(1)*x^2+p(2)*x + p(3);  
    p_new= feval(f,45.8741) + MS_diff_cat_age(roi, 1);
    try
        if max(roots(p))<150, ageing_metric(r, 3)=max(roots(p)); elseif min(roots(p))>30 ageing_metric(r, 3)=min(roots(p)); end; 
    end
    
    ageing_metric(r, 3)=max(roots(p));
    ageing_m_diff(r,1)=ageing_metric(r, 1)-25.514;
    ageing_m_diff(r,2)=ageing_metric(r, 2)-35.4157;
    ageing_m_diff(r,3)=ageing_metric(r, 3)-45.8741;
    
  %f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
  %
  %figure; fplot(f, [17 70], 'LineWidth', 3); %set(gca,'YTickLabel',[]); xlim([18 68]); set(gca,'xtick',[]); hold off; %
end

%https://uk.mathworks.com/help/matlab/ref/roots.html %https://uk.mathworks.com/help/matlab/ref/fzero.html
for r=1:length(id(FDR<0.01))
    roi=id(FDR<0.01); %243% 1:360;% 
    roi=roi(r);
    
    p=p_all(:,roi); f = @(x) p(1)*x^2+p(2)*x + p(3);  
    p_new= feval(f,25.2145) + MS_diff_cat_HC(roi, 1);
    p(3)=p(3)-p_new;
    try
        if max(roots(p))<150, ageing_metric(r, 1)=max(roots(p)); elseif min(roots(p))>30 ageing_metric(r, 1)=min(roots(p)); else ageing_metric(r,3)=0; end; 
    end
    
    p=p_all(:,roi); f = @(x) p(1)*x^2+p(2)*x + p(3);  
    p_new= feval(f,35.6275) + MS_diff_cat_HC(roi, 1);
    p(3)=p(3)-p_new;
    try
        if max(roots(p))<150, ageing_metric(r, 2)=max(roots(p)); elseif min(roots(p))>30 ageing_metric(r, 2)=min(roots(p)); else ageing_metric(r,3)=0; end; 
    end
    
    p=p_all(:,roi); f = @(x) p(1)*x^2+p(2)*x + p(3);  
    p_new= feval(f,47.3538) + MS_diff_cat_HC(roi, 1);
    p(3)=p(3)-p_new;
    try
        if max(roots(p))<150, ageing_metric(r, 3)=max(roots(p)); elseif min(roots(p))>30 ageing_metric(r,3)=min(roots(p)); else ageing_metric(r,3)=0; end; 
    end
    
    ageing_m_diff(r,1)=ageing_metric(r, 1)-25.2145;
    ageing_m_diff(r,2)=ageing_metric(r, 2)-35.6275; 
    ageing_m_diff(r,3)=ageing_metric(r, 3)-47.3538;
    
  %f = @(x) p(1)*x^2+p(2)*x + p(3); %f = @(x) p(1)*x^3+p(2)*x^2 + p(3)*x+ p(4); %
  %
  %figure; fplot(f, [17 70], 'LineWidth', 3); %set(gca,'YTickLabel',[]); xlim([18 68]); set(gca,'xtick',[]); hold off; %
end


clear ageing_m*


clear cluster_val years60 years42
for r=1:length(id(FDR<0.01))
    roi=id(FDR<0.01); %243%1:360
    roi=roi(r);
    p=p_all(:,roi);
    f = @(x) p(1)*x^2+p(2)*x + p(3);
    df = @(x) p(1)*2*x+p(2);
    cluster_val(r, 1)=feval(df, 42.8025); %4th bin
    cluster_val(r, 2)=feval(df, 49.1819); %5th bin
    cluster_val(r, 3)=feval(df, 60); %4th bin
    %MSdiff = m* years >>> years=MSdiff/m
    years1(r,1)=MS_diff_cat_age(roi,1)/cluster_val(r, 3);
    years2(r,1)=MS_diff_cat_age(roi,2)/cluster_val(r, 3);
    years3(r,1)=MS_diff_cat_age(roi,3)/cluster_val(r, 3);
    years4(r,1)=MS_diff_cat_age(roi,4)/cluster_val(r, 3);
    years5(r,1)=MS_diff_cat_age(roi,5)/cluster_val(r, 3);
end

%% updated trajectories
dfx0((tstat4<0 & (T==3) & FDR1<0.01) | (tstat4>0 & (T==4) & FDR1<0.01)& adjr2>0) %c 
label_names_all((tstat4>0 & (T==3) & FDR1<0.01) | (tstat4<0 & (T==4) & FDR1<0.01)& adjr2>0) %c 


%group 1
consistent1=label_names_all(((tstat1>0 & (T==1 |T==3 & dfx0>30) & FDR1<0.01) | (tstat1<0 & (T==2 |T==4 & dfx0>30) & FDR1<0.01)& adjr2>0)... %c
    | ((tstat1>0 & (T==3 & dfx0>30) & FDR1<0.01) | (tstat1<0 & (T==4 & dfx0>30) & FDR1<0.01)& adjr2>0))
inconsistent1=label_names_all(((tstat1<0 & (T==1 |T==3 & dfx0>30) & FDR1<0.01) | (tstat1>0 & (T==2 |T==4 & dfx0>30) & FDR1<0.01)& adjr2>0)... %c
    | ((tstat1<0 & (T==3 & dfx0>30) & FDR1<0.01) | (tstat1>0 & (T==4 & dfx0>30) & FDR1<0.01)& adjr2>0))
%group 2
consistent2=label_names_all(((tstat2>0 & T==1 & FDR2<0.01) | (tstat2<0 & T==2 & FDR2<0.01)& adjr2>0) ...
    |((tstat2>0 & T==4 & dfx0<30 & FDR2<0.01) | (tstat2<0 & T==3 & dfx0<30 & FDR2<0.01)& adjr2>0)...
    |((tstat2<0 & T==4 & dfx0>35 & FDR2<0.01) | (tstat2>0 & T==3 & dfx0>35 & FDR2<0.01)& adjr2>0)) %c
inconsistent2=label_names_all(((tstat2<0 & T==1 & FDR2<0.01) | (tstat2 & T==2 & FDR2<0.01)& adjr2>0) ...
    |((tstat2<0 & T==4 & dfx0<30 & FDR2<0.01) | (tstat2>0 & T==3 & dfx0<30 & FDR2<0.01)& adjr2>0)...
    |((tstat2>0 & T==4 & dfx0>35 & FDR2<0.01) | (tstat2<0 & T==3 & dfx0>35 & FDR2<0.01)& adjr2>0)) %c
%group 3
consistent3=label_names_all(((tstat3>0 & T==1 & FDR3<0.01) | (tstat3<0 & T==2 & FDR3<0.01)& adjr2>0) ...
    |((tstat3>0 & T==4 & dfx0<35 & FDR3<0.01) | (tstat3<0 & T==3 & dfx0<35 & FDR3<0.01)& adjr2>0)...
    |((tstat3<0 & T==4 & dfx0>40 & FDR3<0.01) | (tstat3>0 & T==3 & dfx0>40 & FDR3<0.01)& adjr2>0)) %c
inconsistent3=label_names_all(((tstat3<0 & T==1 & FDR3<0.01) | (tstat3 & T==2 & FDR3<0.01)& adjr2>0) ...
    |((tstat3<0 & T==4 & dfx0<35 & FDR3<0.01) | (tstat3>0 & T==3 & dfx0<35 & FDR3<0.01)& adjr2>0)...
    |((tstat3>0 & T==4 & dfx0>40 & FDR3<0.01) | (tstat3<0 & T==3 & dfx0>40 & FDR3<0.01)& adjr2>0)) %c
%group 4
consistent4=label_names_all(((tstat4>0 & T==1 & FDR4<0.01) | (tstat4<0 & T==2 & FDR4<0.01)& adjr2>0) ...
    |((tstat4>0 & T==4 & dfx0<40 & FDR4<0.01) | (tstat4<0 & T==3 & dfx0<40 & FDR4<0.01)& adjr2>0)...
    |((tstat4<0 & T==4 & dfx0>45 & FDR4<0.01) | (tstat4>0 & T==3 & dfx0>45 & FDR4<0.01)& adjr2>0)) %c
inconsistent4=label_names_all(((tstat4<0 & T==1 & FDR4<0.01) | (tstat4>0 & T==2 & FDR4<0.01)& adjr2>0) ...
    |((tstat4<0 & T==4 & dfx0<40 & FDR4<0.01) | (tstat4>0 & T==3 & dfx0<40 & FDR4<0.01)& adjr2>0)...
    |((tstat4>0 & T==4 & dfx0>45 & FDR4<0.01) | (tstat4<0 & T==3 & dfx0>45 & FDR4<0.01)& adjr2>0)) %c
%group 5
consistent5=label_names_all(((tstat5>0 & T==1 & FDR5<0.01) | (tstat5<0 & T==2 & FDR5<0.01)& adjr2>0) ...
    |((tstat5>0 & T==4 & dfx0<45 & FDR5<0.01) | (tstat5<0 & T==3 & dfx0<45 & FDR5<0.01)& adjr2>0)...
    |((tstat5<0 & T==4 & dfx0>60 & FDR5<0.01) | (tstat5>0 & T==3 & dfx0>60 & FDR5<0.01)& adjr2>0)) %c
inconsistent5=label_names_all(((tstat5<0 & T==1 & FDR5<0.01) | (tstat5>0 & T==2 & FDR5<0.01)& adjr2>0) ...
    |((tstat5<0 & T==4 & dfx0<45 & FDR5<0.01) | (tstat5>0 & T==3 & dfx0<45 & FDR5<0.01)& adjr2>0)...
    |((tstat5>0 & T==4 & dfx0>60 & FDR5<0.01) | (tstat5<0 & T==3 & dfx0>60 & FDR5<0.01)& adjr2>0)) %c


%group 1
consistent1=(((tstat1>0 & (T==1 |T==3 & dfx0>30) & FDR1<0.01) | (tstat1<0 & (T==2 |T==4 & dfx0>30) & FDR1<0.01)& adjr2>0)... %c
    | ((tstat1>0 & (T==3 & dfx0>30) & FDR1<0.01) | (tstat1<0 & (T==4 & dfx0>30) & FDR1<0.01)& adjr2>0))
inconsistent1=(((tstat1<0 & (T==1 |T==3 & dfx0>30) & FDR1<0.01) | (tstat1>0 & (T==2 |T==4 & dfx0>30) & FDR1<0.01)& adjr2>0)... %c
    | ((tstat1<0 & (T==3 & dfx0>30) & FDR1<0.01) | (tstat1>0 & (T==4 & dfx0>30) & FDR1<0.01)& adjr2>0))
%group 2
consistent2=(((tstat2>0 & T==1 & FDR2<0.01) | (tstat2<0 & T==2 & FDR2<0.01)& adjr2>0) ...
    |((tstat2>0 & T==4 & dfx0<30 & FDR2<0.01) | (tstat2<0 & T==3 & dfx0<30 & FDR2<0.01)& adjr2>0)...
    |((tstat2<0 & T==4 & dfx0>35 & FDR2<0.01) | (tstat2>0 & T==3 & dfx0>35 & FDR2<0.01)& adjr2>0)) %c
inconsistent2=(((tstat2<0 & T==1 & FDR2<0.01) | (tstat2 & T==2 & FDR2<0.01)& adjr2>0) ...
    |((tstat2<0 & T==4 & dfx0<30 & FDR2<0.01) | (tstat2>0 & T==3 & dfx0<30 & FDR2<0.01)& adjr2>0)...
    |((tstat2>0 & T==4 & dfx0>35 & FDR2<0.01) | (tstat2<0 & T==3 & dfx0>35 & FDR2<0.01)& adjr2>0)) %c
%group 3
consistent3=(((tstat3>0 & T==1 & FDR3<0.01) | (tstat3<0 & T==2 & FDR3<0.01)& adjr2>0) ...
    |((tstat3>0 & T==4 & dfx0<35 & FDR3<0.01) | (tstat3<0 & T==3 & dfx0<35 & FDR3<0.01)& adjr2>0)...
    |((tstat3<0 & T==4 & dfx0>40 & FDR3<0.01) | (tstat3>0 & T==3 & dfx0>40 & FDR3<0.01)& adjr2>0)) %c
inconsistent3=(((tstat3<0 & T==1 & FDR3<0.01) | (tstat3 & T==2 & FDR3<0.01)& adjr2>0) ...
    |((tstat3<0 & T==4 & dfx0<35 & FDR3<0.01) | (tstat3>0 & T==3 & dfx0<35 & FDR3<0.01)& adjr2>0)...
    |((tstat3>0 & T==4 & dfx0>40 & FDR3<0.01) | (tstat3<0 & T==3 & dfx0>40 & FDR3<0.01)& adjr2>0)) %c
%group 4
consistent4=(((tstat4>0 & T==1 & FDR4<0.01) | (tstat4<0 & T==2 & FDR4<0.01)& adjr2>0) ...
    |((tstat4>0 & T==4 & dfx0<40 & FDR4<0.01) | (tstat4<0 & T==3 & dfx0<40 & FDR4<0.01)& adjr2>0)...
    |((tstat4<0 & T==4 & dfx0>45 & FDR4<0.01) | (tstat4>0 & T==3 & dfx0>45 & FDR4<0.01)& adjr2>0)) %c
inconsistent4=(((tstat4<0 & T==1 & FDR4<0.01) | (tstat4>0 & T==2 & FDR4<0.01)& adjr2>0) ...
    |((tstat4<0 & T==4 & dfx0<40 & FDR4<0.01) | (tstat4>0 & T==3 & dfx0<40 & FDR4<0.01)& adjr2>0)...
    |((tstat4>0 & T==4 & dfx0>45 & FDR4<0.01) | (tstat4<0 & T==3 & dfx0>45 & FDR4<0.01)& adjr2>0)) %c
%group 5
consistent5=(((tstat5>0 & T==1 & FDR5<0.01) | (tstat5<0 & T==2 & FDR5<0.01)& adjr2>0) ...
    |((tstat5>0 & T==4 & dfx0<45 & FDR5<0.01) | (tstat5<0 & T==3 & dfx0<45 & FDR5<0.01)& adjr2>0)...
    |((tstat5<0 & T==4 & dfx0>60 & FDR5<0.01) | (tstat5>0 & T==3 & dfx0>60 & FDR5<0.01)& adjr2>0)) %c
inconsistent5=(((tstat5<0 & T==1 & FDR5<0.01) | (tstat5>0 & T==2 & FDR5<0.01)& adjr2>0) ...
    |((tstat5<0 & T==4 & dfx0<45 & FDR5<0.01) | (tstat5>0 & T==3 & dfx0<45 & FDR5<0.01)& adjr2>0)...
    |((tstat5>0 & T==4 & dfx0>60 & FDR5<0.01) | (tstat5<0 & T==3 & dfx0>60 & FDR5<0.01)& adjr2>0)) %c

tmpcon=sum(sum(consistent1)+sum(consistent2)+sum(consistent3)+sum(consistent4)+sum(consistent5))
tmpinc=sum(sum(inconsistent1)+sum(inconsistent2)+sum(inconsistent3)+sum(inconsistent4)+sum(inconsistent5))
tmpcon/(tmpcon+tmpinc)
