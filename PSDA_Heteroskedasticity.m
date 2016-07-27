%% Dr - Dependent, Sa - Independent
clr
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\';
Sa = importdata(strcat(file_path,'Verify_Sa.txt'));
%Sa = Sa(40:80);
Dr = importdata(strcat(file_path,'Verify_Dr.txt'));
lin_reg = polyfit(log(Sa),log(Dr),1);
homo_mean = lin_reg;
res = log(Dr)-(lin_reg(1)*log(Sa)+lin_reg(2));
subplot(3,2,1)
scatter(log(Sa),log(Dr))
hold on
plot(log(Sa),lin_reg(1)*log(Sa)+lin_reg(2),'linewidth',2)
xlabel('log(Sa)')
ylabel('log(Drift)')
title('Homoskedastic linear regression line and data points')
hold off
grid on
subplot(3,2,2)
scatter(log(Sa),res)
hold on
plot(log(Sa),zeros(max(size(Sa))),'linewidth',2)
hold off
xlabel('log(Sa)')
ylabel('Residuals')
title('Residuals from linear regression under the homoskedastic assumption')
grid on
vari = std(res)^2*ones(max(size(res)),1);
homo_std = sqrt(vari);
Zscore = (log(Dr)-(lin_reg(1)*log(Sa)+lin_reg(2))).^2./vari;
deviance2 = -2*log(-0.5*(sum(log(vari))+sum(Zscore)));
deviance1 = 0;
iter = 1;
while abs(deviance2-deviance1)>=0.01
    iter = iter+1;
param = glmfit([log(Sa)],res.^2,'gamma','log');
vari = glmval(param,[log(Sa)],'log');
lin_reg = lscov([log(Sa) ones(size(Sa))],log(Dr),1./vari);
Zscore = (log(Dr)-(lin_reg(1)*log(Sa)+lin_reg(2))).^2./vari;
deviance1 = deviance2;
deviance2 = -2*log(-0.5*(sum(log(vari))+sum(Zscore)));
end
subplot(3,2,3)
scatter(Dr,sqrt(vari))
hold on
scatter(Dr,homo_std)
hold off
legend('Heteroskedastic','Homoskedastic')
xlabel('Roof drift')
ylabel('Standard deviation')
title('Standard deviation variation with roof drift')
grid on
subplot(3,2,4)
scatter(log(Sa),log(Dr))
hold on
hetero = (lin_reg(1)*log(Sa)+lin_reg(2));
homo = (homo_mean(1)*log(Sa)+homo_mean(2));
plot(log(Sa),(homo),log(Sa),(hetero),'linewidth',1)
hold off
legend('Data','Homoskedastic','Heteroskedastic')
xlabel('Log(Sa)')
ylabel('Log(Roof drift)')
grid on
subplot(3,2,5)
acc = 0.01:0.01:1.4;
plot((acc),sqrt(glmval(param,[log(acc)'],'log')),'o')
xlabel('Sa')
ylabel('Standard deviation')
hold on
scatter(Sa,homo_std)
legend('Heteroskedastic','Homoskedastic')
grid on
hold off
%PSDA
acc = 0.01:0.01:1.4;
std_dev_hetero = sqrt(glmval(param,[log(acc)'],'log'));
std_dev_homo = homo_std(1)*ones(max(size(acc)),1);
sa1_temp = log(importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(sa1_temp(:,1),sa1_temp(:,2),log(acc),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
%dsa1 = abs([dsa1 0]);
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_hetero = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(lin_reg(1)*log(acc)+lin_reg(2)))./std_dev_hetero',0,1);
prob_homo = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(homo_mean(1)*log(acc)+homo_mean(2)))./std_dev_homo',0,1);
AFE_hetero(ii) = sum(prob_hetero.*dsa1)*0.01;
AFE_homo(ii) = sum(prob_homo.*dsa1)*0.01;
end
exact = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact.txt');
subplot(3,2,6)
semilogy(dr_range,-log(1-AFE_hetero)/50,dr_range,-log(1-AFE_homo)/50,'linewidth',1.5)
% hold on
% semilogy(dr_range,AFE_hetero,['.','red'],dr_range,AFE_homo,['.','blue'])
% exact_std = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact_std_drifts.txt');
% figure
% scatter(exact_std(:,1),exact_std(:,2))
% hold on
% scatter(Dr,sqrt(vari))
%% Dr - Independent, Sa - Dependent
clr
Sa = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Sa.txt');
%Sa = Sa(40:80);
Dr = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Drift.txt');
%Dr = Dr(40:80);
lin_reg = polyfit(log(Dr),log(Sa),1);
homo_mean = lin_reg;
res = log(Sa)-(lin_reg(1)*log(Dr)+lin_reg(2));
subplot(3,2,1)
scatter(log(Dr),log(Sa))
hold on
plot(log(Dr),lin_reg(1)*log(Dr)+lin_reg(2),'linewidth',1.5)
ylabel('log(Sa)')
xlabel('log(Drift)')
title('Homoskedastic linear regression line and data points')
hold off
grid on
subplot(3,2,2)
scatter(log(Dr),res)
hold on
plot(log(Dr),zeros(max(size(Dr))),'linewidth',1.5)
hold off
ylabel('log(Sa)')
xlabel('Residuals')
title('Residuals from linear regression under the homoskedastic assumption')
grid on
vari = std(res)^2*ones(max(size(res)),1);
homo_std = sqrt(vari);
Zscore = (log(Sa)-(lin_reg(1)*log(Dr)+lin_reg(2))).^2./vari;
deviance2 = -2*log(-0.5*(sum(log(vari))+sum(Zscore)));
deviance1 = 0;
iter = 0;
while abs(deviance2-deviance1)>=0.01
iter = iter+1;
param = glmfit([log(Dr) log(Dr).^2],res.^2,'gamma','log');
vari = glmval(param,[log(Dr) log(Dr).^2],'log');
lin_reg = lscov([log(Dr) ones(size(Dr))],log(Sa),1./vari);
Zscore = (log(Sa)-(lin_reg(1)*log(Dr)+lin_reg(2))).^2./vari;
deviance1 = deviance2;
deviance2 = -2*log(-0.5*(sum(log(vari))+sum(Zscore)));
end
subplot(3,2,3)
scatter(Sa,sqrt(vari))
hold on
scatter(Sa,homo_std)
hold off
legend('Heteroskedastic','Homoskedastic')
xlabel('Sa')
ylabel('Standard deviation')
title('Standard deviation variation with Sa')
grid on
subplot(3,2,4)
scatter(log(Dr),log(Sa))
hold on
hetero = (lin_reg(1)*log(Dr)+lin_reg(2));
homo = (homo_mean(1)*log(Dr)+homo_mean(2));
plot(log(Dr),(homo),log(Dr),(hetero),'linewidth',1.5)
hold off
legend('Data','Homoskedastic','Heteroskedastic')
ylabel('log(Sa)')
xlabel('log(Roof drift)')
title('Linear regression in log-log space')
grid on
subplot(3,2,5)
plot((Dr),sqrt(glmval(param,[log(Dr) log(Dr).^2],'log')),'o')
xlabel('Roof Drift')
ylabel('Standard deviation')
title('Standard deviation variation with roof drift');
hold on
scatter(Dr,homo_std)
legend('Heteroskedastic','Homoskedastic')
grid on
hold off
%PSDA
acc = 0.01:0.01:2;
std_dev_hetero = sqrt(glmval(param,log(acc),'log'));
std_dev_homo = homo_std(1)*ones(max(size(acc)),1);
sa1_temp = log(importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(sa1_temp(:,1),sa1_temp(:,2),log(acc),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
%dsa1 = abs([dsa1 0]);
dr_range = 0.001:0.001:0.15;
for ii = 1:max(size(dr_range))
prob_hetero = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(lin_reg(1)*log(acc)+lin_reg(2)))./std_dev_hetero',0,1);
prob_homo = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(homo_mean(1)*log(acc)+homo_mean(2)))./std_dev_homo',0,1);
AFE_hetero(ii) = sum(prob_hetero.*dsa1)*0.01;
AFE_homo(ii) = sum(prob_homo.*dsa1)*0.01;
end
% exact = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact.txt');
% figure
% semilogy(dr_range,AFE_hetero,'red',dr_range,AFE_homo,'blue',exact(:,1),exact(:,2),'green')
% hold on
% semilogy(dr_range,AFE_hetero,['.','red'],dr_range,AFE_homo,['.','blue'],exact(:,1),exact(:,2),['.','green'])
exact_std = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact_std_drifts.txt');
subplot(3,2,6)
plot(exact_std(:,1),exact_std(:,2),'red','linewidth',1.5)
title('Comparison of standard deviations')
hold on
plot(sort(Dr),sort(sqrt(vari)),'linewidth',1.5)
xlabel('Roof drift')
ylabel('Standard deviation')
grid on
legend('Exact (~1000 simulations)','Fisher scoring algorithm (88 simulations)');



acc = 0.01:0.01:2;
dr_range = 0.001:0.001:0.20;
std_dev_hetero = sqrt(glmval(param,log(dr_range),'log'));
std_dev_homo = homo_std(1)*ones(max(size(acc)),1);
sa1_temp = log(importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(sa1_temp(:,1),sa1_temp(:,2),log(acc),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
%dsa1 = abs([dsa1 0]);
for ii = 1:max(size(dr_range))
mean_hetero(ii) = (lin_reg(1)*log(dr_range(ii))+lin_reg(2));%1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(homo_mean(1)*log(acc)+homo_mean(2)))./std_dev_hetero(ii),0,1);
mean_homo = (homo_mean(1)*log(dr_range(ii))+homo_mean(2));%1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(homo_mean(1)*log(acc)+homo_mean(2)))./std_dev_homo',0,1);
prob_hetero = normcdf(log(acc),mean_hetero(ii),std_dev_hetero(ii));
prob_homo = normcdf(log(acc),mean_homo,std_dev_homo(ii));
AFE_hetero(ii) = sum(prob_hetero.*dsa1)*0.01;
AFE_homo(ii) = sum(prob_homo.*dsa1)*0.01;
end
mean_homo = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact_mean_drifts.txt');
%plot(dr_range,AFE_hetero)
figure
semilogy(dr_range,-log(1-AFE_hetero)/50,'red',dr_range,-log(1-AFE_homo)/50)
xlabel('Roof drift')
ylabel('AFE')
legend('Heteroskedastic','Homoskedastic')
%% Dr - Dependent, Sa - Independent, New regression
clr
Sa = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Sa.txt');
%Sa = Sa(1:40);
Dr = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Drift.txt');
%Dr = Dr(1:40);
lin_reg = regress(log(Dr),[ones(max(size(Sa)),1) Sa log(Sa)]);
homo_mean = lin_reg;
res = log(Dr)-(lin_reg(1)+lin_reg(2)*Sa+lin_reg(3)*log(Sa));
subplot(3,2,1)
scatter(log(Sa),log(Dr))
hold on
plot(sort(log(Sa)),sort(lin_reg(1)+lin_reg(2)*Sa+lin_reg(3)*log(Sa)),'linewidth',1.5)
xlabel('log(Sa)')
ylabel('log(Drift)')
title('Homoskedastic linear regression line and data points')
hold off
grid on
subplot(3,2,2)
scatter(log(Sa),res)
hold on
plot(log(Sa),zeros(max(size(Sa))),'linewidth',1.5)
hold off
xlabel('log(Sa)')
ylabel('Residuals')
title('Residuals from linear regression under the homoskedastic assumption')
grid on
vari = std(res)^2*ones(max(size(res)),1);
homo_std = sqrt(vari);
Zscore = (log(Dr)-(lin_reg(1)+lin_reg(2)*Sa+lin_reg(3)*log(Sa))).^2./vari;
deviance2 = -2*log(-0.5*(sum(log(vari))+sum(Zscore)));
deviance1 = 0;
iter = 1;
while abs(deviance2-deviance1)>=0.01
param = glmfit(log(Sa),res.^2,'gamma','log');
vari = glmval(param,log(Sa),'log');
lin_reg = lscov([ones(max(size(Dr)),1) Sa log(Sa)],log(Dr),1./vari);
Zscore = (log(Dr)-(lin_reg(1)+lin_reg(2)*Sa+lin_reg(3)*log(Sa))).^2./vari;
deviance1 = deviance2;
deviance2 = -2*log(-0.5*(sum(log(vari))+sum(Zscore)));
iter = iter+1;
end
subplot(3,2,3)
scatter(Dr,sqrt(vari))
hold on
scatter(Dr,homo_std)
hold off
legend('Heteroskedastic','Homoskedastic')
xlabel('Roof drift')
ylabel('Standard deviation')
title('Standard deviation variation with roof drift')
grid on
subplot(3,2,4)
scatter(log(Sa),log(Dr))
hold on
hetero = (lin_reg(1)+lin_reg(2)*Sa+lin_reg(3)*log(Sa));
homo = (homo_mean(3)*log(Sa)+homo_mean(1)+homo_mean(2)*Sa);
plot(sort(log(Sa)),sort((homo)),sort(log(Sa)),sort((hetero)),'linewidth',1.5)
hold off
legend('Data','Homoskedastic','Heteroskedastic')
xlabel('Log(Sa)')
ylabel('Log(Roof drift)')
grid on
subplot(3,2,5)
plot((Sa),sqrt(glmval(param,log(Sa),'log')),'o')
xlabel('Sa')
ylabel('Standard deviation')
hold on
scatter(Sa,homo_std)
legend('Heteroskedastic','Homoskedastic')
grid on
hold off
% PSDA
acc = 0.01:0.01:2;
std_dev_hetero = sqrt(glmval(param,log(acc),'log'));
std_dev_homo = homo_std(1)*ones(max(size(acc)),1);
sa1_temp = log(importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(sa1_temp(:,1),sa1_temp(:,2),log(acc),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
%dsa1 = abs([dsa1 0]);
dr_range = 0.001:0.001:0.15;
for ii = 1:max(size(dr_range))
prob_hetero = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(lin_reg(3)*log(acc)+lin_reg(1)+lin_reg(2)*acc))./std_dev_hetero',0,1);
prob_homo = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(homo_mean(3)*log(acc)+homo_mean(1)+homo_mean(2)*acc))./std_dev_homo',0,1);
AFE_hetero(ii) = sum(prob_hetero.*dsa1)*0.01;
AFE_homo(ii) = sum(prob_homo.*dsa1)*0.01;
end
exact = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact.txt');
subplot(3,2,6)
semilogy(dr_range,-log(1-AFE_hetero)/50,dr_range,-log(1-AFE_homo)/50,'linewidth',1.5)
% hold on
% semilogy(dr_range,AFE_hetero,['.','red'],dr_range,AFE_homo,['.','blue'])
% exact_std = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact_std_drifts.txt');
% figure
% scatter(exact_std(:,1),exact_std(:,2))
% hold on
% scatter(Dr,sqrt(vari))
figure()
ex = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Exact_mean_drifts.txt');
plot(log(ex(:,1)),ex(:,2),sort(hetero),sort(log(Sa)))
%% PSDA
clr
acc = 0.01:0.01:2;
std_dev_220 = sqrt(exp(-2.5133+0.3935*log(acc)));
std_dev_88 = sqrt(exp(-2.3095+0.6428*log(acc)));
std_dev_40 = sqrt(exp(-2.3243+0.2851*log(acc)));
sa1_temp = log(importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(sa1_temp(:,1),sa1_temp(:,2),log(acc),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
%dsa1 = abs([dsa1 0]);
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_220 = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(0.7848*log(acc)-3.4546))./std_dev_220,0,1);
AFE_220(ii) = sum(prob_220.*dsa1)*0.01;

prob_88 = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(0.7421*log(acc)-3.5109))./std_dev_88,0,1);
AFE_88(ii) = sum(prob_88.*dsa1)*0.01;

prob_40 = 1-normcdf((log(dr_range(ii))*ones(max(size(acc)),1)'-(0.7781*log(acc)-3.4712))./std_dev_40,0,1);
AFE_40(ii) = sum(prob_40.*dsa1)*0.01;
end
AFE_220 = -log(1-AFE_220)/50;
AFE_88 = -log(1-AFE_88)/50;
AFE_40 = -log(1-AFE_40)/50;
semilogy(dr_range,AFE_220,'blue',dr_range,AFE_88,'red',dr_range,AFE_40,'green','linewidth',1.5)