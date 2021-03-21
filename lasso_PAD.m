
X_train=MS(group>=4,:); % use rockland sample to train the cross-validated lasso regression 
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

X_HC=MS(group==1,:); % test on HC
pred_age_HC=X_HC*coef + coef0;
corr(pred_age_HC, age_HC)
mean(pred_age_HC - age_HC)
% sqrt(mean((pred_age_HC - age_HC).^2))%boxplot(abs(pred_age_HC - age_HC),'PlotStyle','compact');ylim([0 40]); 

X_SDI=MS(group==2,:); % test on SUD participants, coded as SDI=stimulant dependent individuals
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
% violin plot the brain age gaps
figure(1); violin((pred_age_HC - age_HC), 'facecolor', [0 0 1]); ylim([-30 30])
figure(2); violin(pred_age_SDI-age_SDI, 'facecolor', [1 0 0]); ylim([-30 30])

figure(3); subplot(2,2,1);violin((pred_age_HC(age_HC<30) - age_HC(age_HC<30)), 'facecolor', [0 0 1]); ylim([-30 40]); set(gca,'xtick',[]);set(gca,'ytick',[])
subplot(2,2,2);violin(pred_age_SDI(age_SDI<30)-age_SDI(age_SDI<30), 'facecolor', [1 0 0]); ylim([-30 40]); set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(2,2,3);violin((pred_age_HC(age_HC<45 & age_HC>40) - age_HC(age_HC<45 & age_HC>40)), 'facecolor', [0 0 1]); ylim([-30 40]); set(gca,'xtick',[]);set(gca,'ytick',[])
subplot(2,2,4);violin(pred_age_SDI(age_SDI<45 & age_SDI>40)-age_SDI(age_SDI<45 & age_SDI>40), 'facecolor', [1 0 0]); ylim([-30 40]); set(gca,'xtick',[]); set(gca,'ytick',[])


% compare brain age gaps between HC and SUD groups
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

