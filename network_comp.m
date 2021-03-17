%Ageing_Rmat=Rmat
%ATX_Rmat=Rmat
%ENDO_Rmat=Rmat
%GSK_Rmat=Rmat
    %GENERATE All_Rmat:
        %All_Rmat=horzcat(Ageing_Rmat, ATX_Rmat, ENDO_Rmat, GSK_Rmat)
label_names_HOA={'Frontal Pole (anterior)';'Frontal Pole (ventral)';'Paracingulate Gyrus (anterior)';'Frontal Medial Cortex ';'Frontal Orbital Cortex/Frontal Pole (medial)';'Frontal Orbital Cortex (Medial)';'Postcentral Gyrus';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Subcallosal Cortex (Subgenual ACC)';'Postcentral Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Posterior Cingulate Cortex';'Anterior Cingulate Cortex';'Middle Postcentral Gyrus';'Precentral/postcentral Gyrus';'Precentral Gyrus (lateral ventral)';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Middle Frontal Gyrus';'Frontal Orbital Cortex';'Frontal Orbital Cortex';'Frontal Orbital Cortex (Medial)';'Postcentral Gyrus';'Heschl''s Gyrus';'Middle Frontal Gyrus';'Superior Parietal Lobule';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Precentral Gyrus (medial)';'Precentral Gyrus';'Superior Frontal Gyrus';'Precentral Gyrus';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus';'Superior Parietal Lobule';'Precuneous Cortex ';'Precuneous ';'Superior Parietal Lobule ';'Precuneous Cortex ';'Precuneous ';'Superior/Middle Frontal Gyrus';'Middle Frontal Gyrus';'Frontal Pole';'Medial Superior Frontal Gyrus';'Middle Frontal Gyrus/ Inferior Frontal Gyrus';'Frontal Pole';'Frontal Pole (dorsal anterior)';'Superior Frontal Gyrus/Frontal Pole';'Frontal Pole (anterior)';'Frontal Pole (ventral)';'Heschl''s Gyrus';'Anterior Cingulate Cortex (ventral, anterior)';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate';'Frontal Orbital Cortex/Frontal Pole (lateral)';'Planum Temporale (STG)';'Superior Temporal Gyrus (anterior)';'Middle Frontal Gyrus';'Insular Cortex (anterior)';'Superior Parietal Lobule';'Insular Cortex (anterior)/Frontal Orbital Cortex';'Posterior Cingulate Cortex';'Paracingulate Gyrus';'Precuneous';'Parahippocampal Gyrus';'Precentral Gyrus';'Temporal/Occipital Fisiform Cortex';'Central Opercular Cortex ';'Insular Cortex (dorsal)';'Central Opercular Cortex ';'Frontal Operculum Cortex ';'Frontal Operculum Cortex ';'Lateral Occipital Cortex, inferior division ';'Parahippocampal Gyrus';'Middle Frontal Gyrus ';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus/Inferior Frontal Gyrus';'Inferior Frontal Gyrus, pars triangularis ';'Inferior Frontal Gyrus/Middle Frontal Gyrus';'Insular Cortex/Parietal Operculum';'Lateral Occipital Cortex, superior division ';'Lateral Occipital Cortex, superior division ';'Supramarginal Gyrus, posterior division ';'Lateral Occipital Cortex ';'Planum Temporale (STG)';'Superior Parietal Lobule';'Superior Parietal Lobule ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex, inferior division ';'Heschl''s Gyrus ';'Lateral Occipital Cortex ';'Insular Cortex (dorsal)';'Occipital Cortex';'Lateral Occipital Cortex ';'Frontal Orbital Cortex (Medial)';'Parietal Opercular Cortex';'Insular/Parietal Opercular Cortex';'Central Opercular Cortex';'Frontal Pole (medial)';'Anterior Cingulate';'Anterior Cingulate Cortex';'Medial Superior Frontal Gyrus';'Anterior Cingulate Cortex (anterior)';'Frontal Pole (lateral)';'Middle Frontal Gyrus';'Planum Temporale (STG)';'Precuneous ';'Parahippocampal Gyrus';'Precentral Gyrus';'Parietal Opercular Cortex';'Angular Gyrus/Supramarginal Gyrus';'Postcentral Gyrus (ventral)';'Supramarginal Gyrus, anterior division ';'Supramarginal Gyrus (anterior)';'Angular Gyrus/Lateral Occipital Cortex';'Lateral Occipital Cortex, superior division ';'Angular Gyrus/Lateral Occipital Cortex';'Parahippocampal/Lingual Gyrus';'Temporal Fusiform Cortex, posterior division ';'Temporal Fusiform Cortex';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Planum Polare';'Insular Cortex (ventral)';'Lateral Occipital Cortex ';'Subcallosal Cortex (Subgenual ACC)';'Insular Cortex';'Insular Cortex (posterior)';'Precuneous ';'Precuneous ';'Parahippocampal Gyrus';'Posterior Cingulate Cortex (ventral)';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Parietal Opercular Cortex';'Posterior Cingulate Cortex';'Paracingulate Gyrus (anterior, ventral)';'Middle Frontal Gyrus ';'Supplementary Motor Cortex';'Superior Frontal Gyrus';'Temporal Pole (aTL)';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus (posterior)';'Middle Temporal Gyrus, posterior division';'Middle Temporal Gyrus (posterior)';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Insular Cortex (ventral)';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus, posterior division ';'Inferior Temporal Gyrus, temporooccipital ';'Inferior Temporal Gyrus, posterior division ';'Inferior Temporal Gyrus, temporooccipital part ';'Temporal Fusiform Cortex, posterior division ';'Temporal Pole';'Inferior Temporal Gyrus, anterior division ';'Middle Temporal Gyrus, temporooccipital part ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, superior division ';'Occipital Cortex';'Posterior Cingulate Cortex';'Occipital Cortex';'Occipital Cortex';'Lateral Occipital Cortex ';'Lateral Occipital Cortex, inferior division ';'Occipital Cortex';'Occipital Cortex';'Lateral Occipital Cortex, inferior division ';'Occipital Cortex';'Lateral Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Lateral Occipital Cortex ';'Lingual Gyrus';'Lingual Gyrus';'Temporal Occipital Fusiform Cortex';'Temporal Occipital Fusiform Cortex ';'Frontal Pole (anterior)';'Frontal Pole (ventral)';'Paracingulate Gyrus (anterior)';'Frontal Medial Cortex ';'Frontal Orbital Cortex/Frontal Pole (medial)';'Frontal Orbital Cortex (Medial)';'Postcentral Gyrus';'Posterior Cingulate Cortex';'Posterior Cingulate Cortex';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Subcallosal Cortex (Subgenual ACC)';'Postcentral Gyrus';'Precuneous Cortex ';'Precuneous Cortex ';'Posterior Cingulate Cortex';'Anterior Cingulate Cortex';'Middle Postcentral Gyrus';'Precentral/postcentral Gyrus';'Precentral Gyrus (lateral ventral)';'Inferior Frontal Gyrus, pars opercularis ';'Inferior Frontal Gyrus, pars opercularis ';'Middle Frontal Gyrus';'Frontal Orbital Cortex';'Frontal Orbital Cortex';'Frontal Orbital Cortex (Medial)';'Postcentral Gyrus';'Heschl''s Gyrus';'Middle Frontal Gyrus';'Superior Parietal Lobule';'Medial Postcentral Gyrus ';'Medial Postcentral Gyrus ';'Precentral Gyrus (medial)';'Precentral Gyrus';'Superior Frontal Gyrus';'Precentral Gyrus';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus';'Superior Parietal Lobule';'Precuneous Cortex ';'Precuneous ';'Superior Parietal Lobule ';'Precuneous Cortex ';'Precuneous ';'Superior/Middle Frontal Gyrus';'Middle Frontal Gyrus';'Frontal Pole';'Medial Superior Frontal Gyrus';'Middle Frontal Gyrus/ Inferior Frontal Gyrus';'Frontal Pole';'Frontal Pole (dorsal anterior)';'Superior Frontal Gyrus/Frontal Pole';'Frontal Pole (anterior)';'Frontal Pole (ventral)';'Heschl''s Gyrus';'Anterior Cingulate Cortex (ventral, anterior)';'Anterior Cingulate Cortex (anterior)';'Anterior Cingulate';'Frontal Orbital Cortex/Frontal Pole (lateral)';'Planum Temporale (STG)';'Superior Temporal Gyrus (anterior)';'Middle Frontal Gyrus';'Insular Cortex (anterior)';'Superior Parietal Lobule';'Insular Cortex (anterior)/Frontal Orbital Cortex';'Posterior Cingulate Cortex';'Paracingulate Gyrus';'Precuneous';'Parahippocampal Gyrus';'Precentral Gyrus';'Temporal/Occipital Fisiform Cortex';'Central Opercular Cortex ';'Insular Cortex (dorsal)';'Central Opercular Cortex ';'Frontal Operculum Cortex ';'Frontal Operculum Cortex ';'Lateral Occipital Cortex, inferior division ';'Parahippocampal Gyrus';'Middle Frontal Gyrus ';'Inferior Frontal Gyrus, pars opercularis ';'Precentral Gyrus/Inferior Frontal Gyrus';'Inferior Frontal Gyrus, pars triangularis ';'Inferior Frontal Gyrus/Middle Frontal Gyrus';'Insular Cortex/Parietal Operculum';'Lateral Occipital Cortex, superior division ';'Lateral Occipital Cortex, superior division ';'Supramarginal Gyrus, posterior division ';'Lateral Occipital Cortex ';'Planum Temporale (STG)';'Superior Parietal Lobule';'Superior Parietal Lobule ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex ';'Lateral Occipital Cortex, inferior division ';'Heschl''s Gyrus ';'Lateral Occipital Cortex ';'Insular Cortex (dorsal)';'Occipital Cortex';'Lateral Occipital Cortex ';'Frontal Orbital Cortex (Medial)';'Parietal Opercular Cortex';'Insular/Parietal Opercular Cortex';'Central Opercular Cortex';'Frontal Pole (medial)';'Anterior Cingulate';'Anterior Cingulate Cortex';'Medial Superior Frontal Gyrus';'Anterior Cingulate Cortex (anterior)';'Frontal Pole (lateral)';'Middle Frontal Gyrus';'Planum Temporale (STG)';'Precuneous ';'Parahippocampal Gyrus';'Precentral Gyrus';'Parietal Opercular Cortex';'Angular Gyrus/Supramarginal Gyrus';'Postcentral Gyrus (ventral)';'Supramarginal Gyrus, anterior division ';'Supramarginal Gyrus (anterior)';'Angular Gyrus/Lateral Occipital Cortex';'Lateral Occipital Cortex, superior division ';'Angular Gyrus/Lateral Occipital Cortex';'Parahippocampal/Lingual Gyrus';'Temporal Fusiform Cortex, posterior division ';'Temporal Fusiform Cortex';'Inferior Temporal Gyrus, temporooccipital part ';'Middle Temporal Gyrus, temporooccipital part ';'Planum Polare';'Insular Cortex (ventral)';'Lateral Occipital Cortex ';'Subcallosal Cortex (Subgenual ACC)';'Insular Cortex';'Insular Cortex (posterior)';'Precuneous ';'Precuneous ';'Parahippocampal Gyrus';'Posterior Cingulate Cortex (ventral)';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Parietal Opercular Cortex';'Posterior Cingulate Cortex';'Paracingulate Gyrus (anterior, ventral)';'Middle Frontal Gyrus ';'Supplementary Motor Cortex';'Superior Frontal Gyrus';'Temporal Pole (aTL)';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus (posterior)';'Middle Temporal Gyrus, posterior division';'Middle Temporal Gyrus (posterior)';'Supramarginal Gyrus (posterior)/Angular Gyrus';'Insular Cortex (ventral)';'Middle Temporal Gyrus (anterior)';'Middle Temporal Gyrus, posterior division ';'Inferior Temporal Gyrus, temporooccipital ';'Inferior Temporal Gyrus, posterior division ';'Inferior Temporal Gyrus, temporooccipital part ';'Temporal Fusiform Cortex, posterior division ';'Temporal Pole';'Inferior Temporal Gyrus, anterior division ';'Middle Temporal Gyrus, temporooccipital part ';'Lateral Occipital Cortex, inferior division ';'Lateral Occipital Cortex, superior division ';'Occipital Cortex';'Posterior Cingulate Cortex';'Occipital Cortex';'Occipital Cortex';'Lateral Occipital Cortex ';'Lateral Occipital Cortex, inferior division ';'Occipital Cortex';'Occipital Cortex';'Lateral Occipital Cortex, inferior division ';'Occipital Cortex';'Lateral Occipital Cortex';'Occipital Cortex';'Occipital Cortex';'Lateral Occipital Cortex ';'Lingual Gyrus';'Lingual Gyrus';'Temporal Occipital Fusiform Cortex';'Temporal Occipital Fusiform Cortex '};
label_names_all=vertcat(llabel_names, rlabel_names) %combine all labels

%% Global MS comparisons
%%% HC vs SDI
mean_MS=mean(MS_nodemean')'; y=mean_MS(group<=2 | (group==4 & age<45 & age>20)) %%% group 1=HC, group2=SDI, group 4=ROCKLAND
X=horzcat(group(group<=2 | (group==4 & age<45 & age>20)), sex(group<=2| (group==4 & age<45 & age>20)), age(group<=2| (group==4 & age<45 & age>20))); 
mdl = fitglm(X,y,'linear')


%% MS holds mean similarity for each person for each ROI
mean_MS=mean(MS(group==1,:)) % group mean similarity of all ROIs!
mean_MS2=mean(MS(group==4,:))
tstat_roi=mean(MS((group==4 & age<60 & age>20),:))
%%%%% HC VS SDI
for roi=1:360
    [h,p, ci, stats]=ttest2(MS(1:148, roi), MS(149:331, roi));
    tstat_roi1(roi)=stats.tstat;
    pstat_roi1(roi)=p;
    mean_MS_HC(roi)=mean(MS(group==1,roi)); % group mean similarity of all ROIs!
    mean_MS_SDI(roi)=mean(MS(group==2,roi)); % group mean similarity of all ROIs!
    [p, diff, ES] = permutationTest(MS(1:148, roi), MS(149:331, roi), 5000);
    p_perm(roi)=p;
    diff_perm(roi)=diff;
    ES_perm(roi)=ES;
end

for roi=1:360
    y=MS(group<=2, roi); X=horzcat(group(group<=2), sex(group<=2), age(group<=2)); mdl = fitglm(X,y,'linear');
    tstat_roi(roi)=mdl.Coefficients.tStat(2);
    pstat_roi(roi)=mdl.Coefficients.pValue(2);
end

corr(tstat_roi', tstat_roi1')


%[b,bint,r,rint,stats] = regress(y,X)

FDR = mafdr(pstat_roi1(1:360)) %pstat generated from t-tests or GLMs with t-tests
label_names_all(FDR<0.01) %pstat generated from t-tests or GLMs with t-tests
label_names_all(p_perm<0.01) %pstat generated from permutation testing/t-tests
label_names_HOA(FDR<0.01)
rlabel_names(FDR<0.05)
id=[1:360]'

%%%% HC VS ROCK
for roi=1:360
    [h,p, ci, stats]=ttest2(MS(group==1, roi), MS((group==4 & age<45 & age>20),roi));
    tstat_roi(roi)=stats.tstat;
    pstat_roi(roi)=p;
end
figure; scatter(mean(MS(group==1, :)), mean(MS(group==2 & age<60 & age>20, :)), 6, 'filled')
corr(horzcat(tstat_roi, tstat_roi1'))

FDR_hc_ixi = mafdr(pstat_roi(1:360))
label_names_all(FDR_hc_ixi<0.001)


%%% SDI VS ROCK

for roi=1:360
    [h,p, ci, stats]=ttest2( MS((group==4 & age<45 & age>20), roi), MS(group==2, roi));
    tstat_roi(roi)=stats.tstat;
    pstat_roi(roi)=p;
end
FDR_sdi_ixi = mafdr(pstat_roi(1:360))
label_names_all(FDR_sdi_ixi<0.001) %& FDR<0.01)
corr(tstat_roi, tstat_roi1')
scatter(tstat_roi', tstat_roi1', 8, 'k', 'filled'); %x=tstat_roi'; X=horzcat(ones(1,360)', x); y=tstat_roi1'; b=x/y

label_names_all((FDR_hc_sdi<0.01) & (FDR_hc_ixi<0.001) & ~(FDR_sdi_ixi<0.001)) %ageing and sdi

label_names_all(~(FDR_hc_sdi<0.01) & (FDR_hc_ixi<0.001) & (FDR_sdi_ixi<0.001)) % uniquely ageing

label_names_all((FDR_hc_sdi<0.01) & ~(FDR_hc_ixi<0.001) & (FDR_sdi_ixi<0.001)) %uniquely sdi


%% INDIVIDUAL datasets:
%ENDO
%clear pstat_roi
for roi=1:360
    [h,p, ci, stats]=ttest2(MS(datasetNr==1 & group==1, roi), MS((datasetNr==1 & group==2), roi));
    tstat_1(roi)=stats.tstat;
    pstat_1(roi)=p;
end


%ATX
for roi=1:360
    [h,p, ci, stats]=ttest2(MS((datasetNr==2 & group==1), roi), MS((datasetNr==2 & group==2), roi));
    tstat_2(roi)=stats.tstat;
    pstat_2(roi)=p;
end


%GSK
for roi=1:360
    [h,p, ci, stats]=ttest2(MS((datasetNr==3 & group==1), roi), MS((datasetNr==3 & group==2), roi));
    tstat_3(roi)=stats.tstat;
    pstat_3(roi)=p;
end
%FDR = mafdr(pstat_roi(1:360))
%label_names_all(FDR<0.05)

%Ageing
for roi=1:360
    [h,p, ci, stats]=ttest2(MS((datasetNr==4 & group==1), roi), MS((datasetNr==4 & group==2), roi));
    tstat_4(roi)=stats.tstat;
    pstat_4(roi)=p;
end

corr(horzcat(tstat_1', tstat_2', tstat_3', tstat_4'))
figure; 
subplot(2,3,1); scatter(tstat_1, tstat_2, 5, 'k', 'filled'); x=tstat_1'; y=tstat_2';
X = [ones(length(x),1) x]; b = X\y; yCalc = X*b; hold on; plot(x,yCalc,'k'); ylabel('ATX'); xlabel('ENDO'); xlim([-4 4]);ylim([-4 4]);
subplot(2,3,2); scatter(tstat_1, tstat_2, 5, 'k', 'filled'); x=tstat_1'; y=tstat_3';
X = [ones(length(x),1) x]; b = X\y; yCalc = X*b; hold on; plot(x,yCalc,'k'); ylabel('GSK'); xlabel('ENDO'); xlim([-4 4]);ylim([-4 4])
subplot(2,3,3); scatter(tstat_1, tstat_2, 5, 'k', 'filled'); x=tstat_1'; y=tstat_4';
X = [ones(length(x),1) x]; b = X\y; yCalc = X*b; hold on; plot(x,yCalc,'k'); ylabel('Ageing'); xlabel('ENDO'); xlim([-4 4]);ylim([-4 4])
subplot(2,3,4); scatter(tstat_1, tstat_2, 5, 'k', 'filled'); x=tstat_2'; y=tstat_3';
X = [ones(length(x),1) x]; b = X\y; yCalc = X*b; hold on; plot(x,yCalc,'k'); ylabel('GSK'); xlabel('ATX'); xlim([-4 4]);ylim([-4 4])
subplot(2,3,5); scatter(tstat_1, tstat_2, 5, 'k', 'filled'); x=tstat_2'; y=tstat_4';
X = [ones(length(x),1) x]; b = X\y; yCalc = X*b; hold on; plot(x,yCalc,'k'); ylabel('Ageing'); xlabel('ATX');xlim([-4 4]);ylim([-4 4])
subplot(2,3,6); scatter(tstat_1, tstat_2, 5, 'k', 'filled'); x=tstat_3'; y=tstat_4';
X = [ones(length(x),1) x]; b = X\y; yCalc = X*b; hold on; plot(x,yCalc,'k'); ylabel('Ageing'); xlabel('GSK'); xlim([-4 4]);ylim([-4 4])

