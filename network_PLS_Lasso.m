clear X Y
X=zscore(MS(group==2,:));
Y(:,1)=zscore(age_SDI); %Y=age(group==2);
permutations=5000;   
allobservations=Y; 

for ncomp=1:8
    
    for n = 1:permutations
    % selecting either next combination, or random permutation
    permutation_index = randperm(length(allobservations));
    % creating random sample based on permutation index
    randomSample = allobservations(permutation_index,:);
    % running the PLS for this permutation
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,randomSample,ncomp);
    Rsq(n) = sum(PCTVAR(2,:));
    Rsq1(n) = sum(PCTVAR(2,:));
    end
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp);
    p(ncomp)=sum(sum(PCTVAR(2,:))<Rsq')/permutations
    p_1(ncomp)=sum(sum(PCTVAR(1,:))<Rsq')/permutations
end
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp);
%PCTVAR  0.0949    0.3619    0.1476    0.0778    0.1011    0.0831    0.0493    0.0283
% P   0.7232    0.4158    0.1454    0.1906    0.2490    0.2168    0.2072    0.2852

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear X Y
X=zscore(MS(group==1,:));
Y(:,1)=zscore(age_HC); %Y=age(group==2);
permutations=5000;   
allobservations=Y; 

for ncomp=1:8
    
    for n = 1:permutations
    % selecting either next combination, or random permutation
    permutation_index = randperm(length(allobservations));
    % creating random sample based on permutation index
    randomSample = allobservations(permutation_index,:);
    % running the PLS for this permutation
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,randomSample,ncomp);
    Rsq(n) = sum(PCTVAR(2,:));
    Rsq1(n) = sum(PCTVAR(2,:));
    end
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp);
    p(ncomp)=sum(sum(PCTVAR(2,:))<Rsq')/permutations
    p_1(ncomp)=sum(sum(PCTVAR(1,:))<Rsq')/permutations
end
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,Y,ncomp);
%PCTVAR 0.5044    0.1228    0.0413    0.0788    0.0936    0.0684    0.0377    0.0231
%P  0.0058    0.0302    0.1978    0.3676    0.3608    0.5042    0.5970    0.6244
corr(XL, tstat_roi1')

%% LASSO

clear X Y
X=zscore(MS(group==2,:));
Y(:,1)=zscore(age_SDI); %Y=age(group==2);
permutations=5000;   
allobservations=Y; 

% basic model: %[B,FitInfo] = lasso(X,Y)
[B_SDI,FitInfo] = lasso(X,Y,'CV',10,'PredictorNames', label_names_all);
idxLambdaMinMSE = FitInfo.IndexMinMSE;
minMSEModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinMSE)~=0)
coef = B(:,idxLambdaMinMSE);
coef0 = FitInfo.Intercept(idxLambdaMinMSE);
yhat = X*coef + coef0;
scatter(Y,yhat); xlim([-3 3]); ylim([-1 1]); hold on; plot(Y,Y); hold off
            %less useful!
            idxLambda1SE = FitInfo.Index1SE;
            sparseModelPredictors = FitInfo.PredictorNames(B(:,idxLambda1SE)~=0)
            lassoPlot(B,FitInfo,'PlotType','CV');
            legend('show') % Show legend

            

            
clear X Y
X=zscore(MS(group==1,:));
Y(:,1)=zscore(age_HC); %Y=age(group==2);
permutations=5000;   
allobservations=Y; 

[B_HC,FitInfo] = lasso(X,Y,'CV',10,'PredictorNames', label_names_all);
idxLambdaMinMSE = FitInfo.IndexMinMSE;
minMSEModelPredictors = FitInfo.PredictorNames(B(:,idxLambdaMinMSE)~=0)
coef = B(:,idxLambdaMinMSE);
coef0 = FitInfo.Intercept(idxLambdaMinMSE);
yhat = X*coef + coef0;
scatter(Y,yhat); xlim([-3 3]); ylim([-1 1]); hold on; plot(Y,Y); hold off

%train the lasso on Rockland sample and then use the model to predict the
%age of the SDI and HC; see how it deviates

X_train=MS(group>=4,:);
Y_train=age(group>=4);
[B_ROCK,FitInfo] = lasso(X_train,Y_train,'CV',10,'PredictorNames', label_names_all);
idxLambdaMinMSE = FitInfo.IndexMinMSE;
minMSEModelPredictors = FitInfo.PredictorNames(B_ROCK(:,idxLambdaMinMSE)~=0)
coef = B_ROCK(:,idxLambdaMinMSE);
coef0 = FitInfo.Intercept(idxLambdaMinMSE);
yhat = X_train*coef + coef0;
scatter(Y_train,yhat);hold on; plot(Y_train,Y_train, 'k'); hold off % xlim([-3 3]); ylim([-1 1]); 
%mean(abs(Y_train-yhat));%boxplot(abs(Y_train-yhat),'PlotStyle','compact');ylim([0 40]); 
corr(Y_train,yhat)

X_HC=MS(group==1,:);
pred_age_HC=X_HC*coef + coef0;
corr(pred_age_HC, age_HC)
mean(pred_age_HC - age_HC)
% sqrt(mean((pred_age_HC - age_HC).^2))%boxplot(abs(pred_age_HC - age_HC),'PlotStyle','compact');ylim([0 40]); 

X_SDI=MS(group==2,:);
pred_age_SDI=X_SDI*coef + coef0;
corr(pred_age_SDI, age_SDI)
mean(pred_age_SDI - age_SDI)
%mean(abs(pred_age_SDI - age_SDI)) %boxplot(abs(pred_age_SDI - age_SDI),'PlotStyle','compact');ylim([0 40]); 

mean(pred_age_SDI(age_SDI<30) - age_SDI(age_SDI<30))
mean(pred_age_SDI(30<age_SDI & age_SDI <35) - age_SDI(30<age_SDI & age_SDI <35))
mean(pred_age_SDI(35<age_SDI & age_SDI <40) - age_SDI(35<age_SDI & age_SDI <40))
mean(pred_age_SDI(40<age_SDI & age_SDI <45) - age_SDI(40<age_SDI & age_SDI <45))
mean(pred_age_SDI(45<age_SDI & age_SDI <60) - age_SDI(45<age_SDI & age_SDI <60))


[h, p, ci, stats] =ttest2((pred_age_HC - age_HC), (pred_age_SDI-age_SDI))
% violin plot
figure(1); violin((pred_age_HC - age_HC), 'facecolor', [0 0 1]); ylim([-30 30])
figure(2); violin(pred_age_SDI-age_SDI, 'facecolor', [1 0 0]); ylim([-30 30])

figure(3); subplot(2,2,1);violin((pred_age_HC(age_HC<30) - age_HC(age_HC<30)), 'facecolor', [0 0 1]); ylim([-30 40]); set(gca,'xtick',[]);set(gca,'ytick',[])
subplot(2,2,2);violin(pred_age_SDI(age_SDI<30)-age_SDI(age_SDI<30), 'facecolor', [1 0 0]); ylim([-30 40]); set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(2,2,3);violin((pred_age_HC(age_HC<45 & age_HC>40) - age_HC(age_HC<45 & age_HC>40)), 'facecolor', [0 0 1]); ylim([-30 40]); set(gca,'xtick',[]);set(gca,'ytick',[])
subplot(2,2,4);violin(pred_age_SDI(age_SDI<45 & age_SDI>40)-age_SDI(age_SDI<45 & age_SDI>40), 'facecolor', [1 0 0]); ylim([-30 40]); set(gca,'xtick',[]); set(gca,'ytick',[])



[h, p, ci, stats] =ttest2((pred_age_HC(age_HC<30) - age_HC(age_HC<30)), (pred_age_SDI(age_SDI<30)-age_SDI(age_SDI<30)))
[h, p, ci, stats] =ttest2((pred_age_HC(30<age_HC & age_HC <35) - age_HC(30<age_HC & age_HC <35)), (pred_age_SDI(30<age_SDI & age_SDI <35)-age_SDI(30<age_SDI & age_SDI <35)))
[h, p, ci, stats] =ttest2((pred_age_HC(35<age_HC & age_HC <40) - age_HC(35<age_HC & age_HC <40)), (pred_age_SDI(35<age_SDI & age_SDI <40)-age_SDI(35<age_SDI & age_SDI <40)))
[h, p, ci, stats] =ttest2((pred_age_HC(40<age_HC & age_HC <45) - age_HC(40<age_HC & age_HC <45)), (pred_age_SDI(40<age_SDI & age_SDI <45)-age_SDI(40<age_SDI & age_SDI <45)))
[h, p, ci, stats] =ttest2((pred_age_HC(45<age_HC & age_HC <60) - age_HC(45<age_HC & age_HC <60)), (pred_age_SDI(45<age_SDI & age_SDI <60)-age_SDI(45<age_SDI & age_SDI <60)))



mean(pred_age_SDI(age_HC<30) - age_SDI(age_HC<30))
mean(pred_age_HC(30<age_HC & age_HC <35) - age_HC(30<age_HC & age_HC <35))
mean(pred_age_HC(35<age_HC & age_HC <40) - age_HC(35<age_HC & age_HC <40))
mean(pred_age_HC(40<age_HC & age_HC <45) - age_HC(40<age_HC & age_HC <45))
mean(pred_age_HC(45<age_HC & age_HC <60) - age_HC(45<age_HC & age_HC <60))


%permutation testing this thing

X_train=MS(group>=4,:);
Y_train=age(group>=4);
8

%% Gaussian regression as a more sensitive age test
X_train=MS(group==4,:);
Y_train=age(group==4);
cvgprMdl = fitrgp(X_train,Y_train ,'KernelFunction','ardsquaredexponential',...
      'FitMethod','sr','PredictMethod','fic','Standardize',1)
L = resubLoss(cvgprMdl)

[ypred, ysd, yint] = resubPredict(cvgprMdl);
plot(ypred); hold on; plot(Y_train)
scatter(ypred, Y_train)
corr(ypred, Y_train)


[ypred, ysd, yint]=predict(cvgprMdl,X_HC);
corr(pred_age_HC, age_HC)
mean(pred_age_HC - age_HC)


X_SDI=MS(group==2,:);
pred_age_SDI=X_SDI*coef + coef0;
corr(pred_age_SDI, age_SDI)
mean(pred_age_SDI - age_SDI)


%% svm
X_train=MS(group==4,:);
Y_train=age(group==4);
Mdl=fitrsvm(X_train,Y_train, 'KernelFunction','gaussian', 'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))
L = resubLoss(Mdl)

ypred = resubPredict(Mdl);
%plot(ypred); hold on; plot(Y_train)
scatter(ypred, Y_train)
corr(ypred, Y_train)

%HC
pred_age_HC=predict(Mdl,X_HC);
corr(pred_age_HC, age_HC)
scatter(pred_age_HC, age_HC)
mean(pred_age_HC - age_HC)

%SDI
pred_age_SDI=predict(Mdl,X_SDI);
corr(pred_age_SDI, age_SDI)
mean(pred_age_SDI - age_SDI)

[h, p, ci, stats] =ttest2((pred_age_HC - age_HC), (pred_age_SDI-age_SDI))

%IXI
pred_age_IXI=predict(Mdl,MS(group==3,:));
corr(pred_age_IXI, age(group==3))
scatter(pred_age_IXI, age(group==3))
mean(pred_age_IXI - age(group==3))


%% alternative - pls

X_train=MS(group==4,:);
Y_train=age(group==4);

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X_train,Y_train,5);

yhat = [ones(size(X_train ,1),1) X_train]*BETA;

scatter(Y_train,yhat);hold on; plot(Y_train,Y_train); hold off % xlim([-3 3]); ylim([-1 1]); 

corr(Y_train,yhat)


X_HC=MS(group==1,:);
pred_age_HC = [ones(size(X_HC,1),1) X_HC]*BETA;

corr(pred_age_HC, age_HC)
mean(pred_age_HC - age_HC)


X_SDI=MS(group==2,:);
pred_age_SDI = [ones(size(X_SDI,1),1) X_SDI]*BETA;
corr(pred_age_SDI, age_SDI)
mean(pred_age_SDI - age_SDI)

mean(pred_age_SDI(age_SDI<30) - age_SDI(age_SDI<30))
mean(pred_age_SDI(30<age_SDI & age_SDI <35) - age_SDI(30<age_SDI & age_SDI <35))
mean(pred_age_SDI(35<age_SDI & age_SDI <40) - age_SDI(35<age_SDI & age_SDI <40))
mean(pred_age_SDI(40<age_SDI & age_SDI <45) - age_SDI(40<age_SDI & age_SDI <45))
mean(pred_age_SDI(45<age_SDI & age_SDI <60) - age_SDI(45<age_SDI & age_SDI <60))
scatter(pred_age_SDI, age_SDI)

[h, p, ci, stats] =ttest2((pred_age_HC - age_HC), (pred_age_SDI-age_SDI))

histogram(corr(MS, age))