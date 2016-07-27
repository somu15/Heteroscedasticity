%% Fisher's scoring algorithm
% Written by Somayajulu on March 18, 2016
% CHANGE THE PATH IF NECESSARY
% Dr - Independent, Sa - Dependent
% This is a more case specific algorithm. Works only when log of dependent 
% variable follows a linear model with normal distribution and log of 
% variance follows a linear model with gamma distribution.
% NOTE: THIS IS NOT THE ORIGINAL FISHER'S SCORING ALGORITHM. IT IS A 
% SIMPLIFIED VERSION.
clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\';
Dr = importdata(strcat(file_path,'Verify_Sa.txt'));
%Sa = Sa(40:80);
Sa = importdata(strcat(file_path,'Verify_Dr.txt'));
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
scatter(log(Sa),res)
hold on
plot(log(Sa),zeros(max(size(Dr))),'linewidth',1.5)
hold off
xlabel('log(Sa)')
ylabel('Residuals')
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
param = glmfit(log(Dr),res.^2,'gamma','log');
vari = glmval(param,log(Dr),'log');
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
plot((Dr),sqrt(vari),'o')
xlabel('Roof Drift')
ylabel('Standard deviation')
title('Standard deviation variation with roof drift');
hold on
scatter(Dr,homo_std)
legend('Heteroskedastic','Homoskedastic')
grid on
hold off
exact_std = importdata(strcat(file_path,'Exact_std_drifts1.txt'));
subplot(3,2,6)
plot(exact_std(:,1),exact_std(:,2),'red','linewidth',1.5)
title('Comparison of standard deviations')
hold on
plot(0.008:0.001:0.08,sqrt(exp(param(1)+param(2)*log(0.008:0.001:0.08))),'linewidth',1.5)
xlabel('Roof drift')
ylabel('Standard deviation')
grid on
legend('Exact (~1000 simulations)','Fisher scoring algorithm (88 simulations)');
%% Fisher's scoring algorithm
% Written by Somayajulu on March 25, 2016
% Modified by Somayajulu on April 1, 2016 
% CHANGE THE PATH IF NECESSARY
% Dr - Dependent, Sa - Independent
% NOTE: THE FIRST PART OF THE CODE IS AITKIN's ORIGINAL FISHER SCORING ALGORITHM
% NOTE: THE SECOND PART OF THE CODE IS GIBBS-METROPOLIS ALGORITHM WITH
% CEPEDA-GAMERMAN CORRECTION FOR DEPENDENT VARIABLE TO ENSURE GOOD
% ACCEPTANCE RATES FOR PROPOSALS (this is Bayesian)
% QUADRATIC VARIANCE FUNCTION 
clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA1\';
Sa = importdata(strcat(file_path,'x.txt'));
% display('Independent is Sa3')
%fp = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\RD\';
% eps = importdata(strcat(fp,'eps_sa1.txt'));
% display('eps is Sa3')
%Sa = Sa(1:41);
Dr = importdata(strcat(file_path,'y.txt'));


% M = importdata(strcat(fp,'Mw.txt'));
% R = importdata(strcat(fp,'Rjb.txt'));
X = [ones(max(size(Sa)),1) log(Sa)];
Z = [ones(max(size(Sa)),1) (Sa) (Sa).^2];
beta_0 = regress(log(Dr),X);
res = log(Dr)-(beta_0(1)+beta_0(2)*log(Sa));
OLS_res = res;
vari = std(res)^2*ones(max(size(Sa)),1);
% str = regstats(res,M,'linear',{'tstat'});
% str.tstat
% str1 = regstats(res,R,'linear',{'tstat'});
% str1.tstat
% str2 = regstats(res,eps,'linear',{'tstat'});
% str2.tstat



OLS_std = sqrt(vari);
gamma_0 = [log(std(res)^2);0;0];
parameters = vertcat(beta_0,gamma_0);
W = diag(1./vari);
A_0 = blkdiag(X'*W*X,0.5*Z'*Z);
temp = X'*(res./vari);
temp1 = 0.5*Z'*(res.^2./vari-1);
d_0 = vertcat(temp,temp1);
count = 1;
delta = 3;
while norm(delta)>=0.01
delta = inv(A_0)*d_0;
parameters = parameters+delta;
res = log(Dr)-(parameters(1)+parameters(2)*log(Sa));
vari = (exp(parameters(3)+parameters(4)*(Sa)+parameters(5)*(Sa).^2));
temp = X'*(res./vari);
temp1 = 0.5*Z'*(res.^2./vari-1);
d_0 = vertcat(temp,temp1);
W = diag(1./vari);
A_0 = blkdiag(X'*W*X,0.5*Z'*Z);
count = count+1;
if count>1000
break
end
end
parameters_Freq = parameters;
parameters_Freq
beta_mean_pri = parameters(1:2);
gamma_mean_pri = 0.0000001*parameters(3:max(size(parameters)));
cov_matrix = inv(A_0);
beta_cov_pri = cov_matrix(1:2,1:2);
gamma_cov_pri = cov_matrix(3:max(size(parameters)),3:max(size(parameters)));
gamma_cov_pri = diag(diag(gamma_cov_pri))+[10 0 0;0 10^10 0;0 0 1];
gamma_old = mvnrnd(gamma_mean_pri,gamma_cov_pri);
vari_old = vari;
count = 0;
max_iter = 10000;
for ii = 1:max_iter
beta_cov_post = inv(inv(beta_cov_pri)+X'*diag(vari_old)*X);
beta_mean_post = beta_cov_post*(inv(beta_cov_pri)*beta_mean_pri+X'*diag(vari_old)*log(Dr));
beta_prop(ii,:) = mvnrnd(beta_mean_post,beta_cov_post);
Trans = log(vari_old)+(log(Dr)-(beta_prop(ii,:)*X')').^2./vari_old-1;
gamma_cov_prop = inv(inv(gamma_cov_pri)+0.5*Z'*Z);
gamma_mean_prop = gamma_cov_prop*(inv(gamma_cov_pri)*gamma_mean_pri+0.5*Z'*Trans);
gamma_prop = mvnrnd(gamma_mean_prop,gamma_cov_prop);
vari_new = exp(gamma_prop(1)+gamma_prop(2)*(Sa)+gamma_prop(3)*(Sa).^2);
accept_ratio = sum(log(normpdf(log(Dr),(beta_prop(ii,:)*X')',sqrt(vari_new)))-(log(normpdf(log(Dr),(beta_prop(ii,:)*X')',sqrt(vari_old)))))+log(mvnpdf(gamma_prop',gamma_mean_pri,gamma_cov_pri))-log(mvnpdf(gamma_old',gamma_mean_pri,gamma_cov_pri));
decision = log(rand);
if decision<accept_ratio
gamma_new(ii,:) = gamma_prop;
vari_old = vari_new;
gamma_old = gamma_prop;
count = count+1;
else
gamma_new(ii,:) = gamma_old;    
end
end
ii = 1:max_iter;
parameters_bayes(1) = mean(beta_prop(:,1));
parameters_bayes(2) = mean(beta_prop(:,2));
parameters_bayes(3) = mean(gamma_new(1000:max_iter,1));
parameters_bayes(4) = mean(gamma_new(1000:max_iter,2));
parameters_bayes(5) = mean(gamma_new(1000:max_iter,3));
cov_matrix_bayes(1:2,1:2) = cov(beta_prop(1000:max_iter,1),beta_prop(1000:max_iter,2));
cov_matrix_bayes(3:4,3:4) = cov(gamma_new(1000:max_iter,1),gamma_new(1000:max_iter,2));
K = cov(gamma_new(1000:max_iter,2),gamma_new(1000:max_iter,3));
cov_matrix_bayes(4,5) = K(1,2);
cov_matrix_bayes(5,4) = K(1,2);
cov_matrix_bayes(5,5) = K(2,2);
K = cov(gamma_new(1000:max_iter,1),gamma_new(1000:max_iter,3));
cov_matrix_bayes(3,5) = K(1,2);
cov_matrix_bayes(5,3) = K(1,2);
parameters_bayes = parameters_bayes';
parameters_bayes
cov_matrix
cov_matrix_bayes
%K = exp(parameters_bayes(3)+parameters_bayes(4)*(Sa)+parameters_bayes(5)*(Sa).^2);
subplot(1,2,2)
hist(Sa)
xlabel('Sa(T1)')
ylabel('Frequency')
subplot(1,2,1)
% subplot(3,2,1)
% plot(10:max_iter,beta_prop(10:max_iter,1))
% grid on
% xlabel('Iteration')
% ylabel('Beta 1')
% subplot(3,2,2)
% plot(10:max_iter,beta_prop(10:max_iter,2))
% grid on
% xlabel('Iteration')
% ylabel('Beta 2')
% subplot(3,2,3)
% plot(1:max_iter,gamma_new(1:max_iter,1))
% grid on
% xlabel('Iteration')
% ylabel('Gamma 1')
% subplot(3,2,4)
% plot(10:max_iter,gamma_new(10:max_iter,2))
% grid on
% xlabel('Iteration')
% ylabel('Gamma 2')
% subplot(3,2,5)
% plot(10:max_iter,gamma_new(10:max_iter,3))
% grid on
% xlabel('Iteration')
% ylabel('Gamma 3')
% subplot(3,2,6)
sa_range = 0.05:0.05:2;
exact = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA1\Exact_stds_2g.txt');
reg_para = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA1\reg.txt');
plot(sa_range,sqrt((exp(parameters(3)+parameters(4)*(sa_range)+parameters(5)*(sa_range).^2))),'red','linewidth',1.5)
hold on
plot(sa_range,sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(sa_range)+parameters_bayes(5)*(sa_range).^2)),sa_range,sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2))),'green','linewidth',1.5)%
hold on
% ex = interp1(exact(:,1),exact(:,2),sa_range);
% scatter(exact(:,1),exact(:,2))
% hold on
scatter(exact(:,1),exact(:,2))
%rBayes = rsquare(exact(:,2),sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(exact(:,1))+parameters_bayes(5)*(exact(:,1)).^2)));
SSE_Bayes = sum((exact(:,2)-sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(exact(:,1))+parameters_bayes(5)*(exact(:,1)).^2))).^2);
%rFreq = rsquare(exact(:,2),sqrt((exp(parameters(3)+parameters(4)*(exact(:,1))+parameters(5)*(exact(:,1)).^2))));
SSE_Freq = sum((exact(:,2)-sqrt((exp(parameters(3)+parameters(4)*(exact(:,1))+parameters(5)*(exact(:,1)).^2)))).^2);
%rFit = rsquare(exact(:,2),sqrt((exp(-3.83+2.91*(exact(:,1))-1.02*(exact(:,1)).^2))));
SSE_Fit = sum((exact(:,2)'-sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2)))).^2);
% rBayes
% rFreq
% rFit
SSE_Bayes
SSE_Freq
SSE_Fit
reg_para
std_freq = sqrt((exp(parameters(3)+parameters(4)*(sa_range)+parameters(5)*(sa_range).^2)));
std_bayes = sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(sa_range)+parameters_bayes(5)*(sa_range).^2));
std_exact = sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2)));
grid on
xlabel('Sa (g)')
ylabel('Standard deviation in roof drift')
legend(strcat('Frequentist SSE = ',num2str(SSE_Freq)),strcat('Bayesian SSE = ',num2str(SSE_Bayes)),strcat('Fit to actual data SSE = ',num2str(SSE_Fit)))
%% Fisher's scoring algorithm
% Written by Somayajulu on March 25, 2016
% CHANGE THE PATH IF NECESSARY
% Dr - Dependent, Sa - Independent
% NOTE: THIS IS THE AITKIN's ORIGINAL SCORING ALGORITHM
% LINEAR VARIANCE FUNCTION 
clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA4\';
Sa = importdata(strcat(file_path,'x.txt'));
%Sa = Sa(1:40);
Dr = importdata(strcat(file_path,'y.txt'));
%Dr = Dr(find(Sa<1.1));
%Sa = Sa(find(Sa<1.1));
X = [ones(max(size(Sa)),1) log(Sa)];
Z = [ones(max(size(Sa)),1) (Sa)];
beta_0 = regress(log(Dr),X);
res = log(Dr)-(beta_0(1)+beta_0(2)*(Sa));
OLS_res = res;
vari = std(res)^2*ones(max(size(Sa)),1);
OLS_std = sqrt(vari);
gamma_0 = [log(std(res)^2);0];%%Log here
parameters = vertcat(beta_0,gamma_0);
W = diag(1./vari);
A_0 = blkdiag(X'*W*X,0.5*Z'*Z);
temp = X'*(res./vari);
temp1 = 0.5*Z'*(res.^2./vari-1);
d_0 = vertcat(temp,temp1);
count = 1;
delta = 3;
while norm(delta)>=0.01
delta = inv(A_0)*d_0;
parameters = parameters+delta;
res = log(Dr)-(parameters(1)+parameters(2)*log(Sa));
vari = (exp(parameters(3)+parameters(4)*(Sa)));%%Exp here
temp = X'*(res./vari);
temp1 = 0.5*Z'*(res.^2./vari-1);
d_0 = vertcat(temp,temp1);
W = diag(1./vari);
A_0 = blkdiag(X'*W*X,0.5*Z'*Z);
count = count+1;
if count>100
break
end
end
parameters
beta_mean_pri = parameters(1:2);
gamma_mean_pri = 0.000000001*parameters(3:max(size(parameters)));
cov_matrix = inv(A_0);
beta_cov_pri = cov_matrix(1:2,1:2);
gamma_cov_pri = cov_matrix(3:max(size(parameters)),3:max(size(parameters)))+[10^10 0;0 1];
Sa_range = 0.01:0.01:1;
%plot((Sa_range),sqrt((exp(parameters(3)+parameters(4)*(Sa_range)))))
gamma_old = mvnrnd(gamma_mean_pri,gamma_cov_pri);
vari_old = vari;
count = 0;
max_iter = 10000;
for ii = 1:max_iter
beta_cov_post = inv(inv(beta_cov_pri)+X'*diag(vari_old)*X);
beta_mean_post = beta_cov_post*(inv(beta_cov_pri)*beta_mean_pri+X'*diag(vari_old)*log(Dr));
beta_prop(ii,:) = mvnrnd(beta_mean_post,beta_cov_post);
Trans = log(vari_old)+(log(Dr)-(beta_prop(ii,:)*X')').^2./vari_old-1;
gamma_cov_prop = inv(inv(gamma_cov_pri)+0.5*Z'*Z);
gamma_mean_prop = gamma_cov_prop*(inv(gamma_cov_pri)*gamma_mean_pri+0.5*Z'*Trans);
gamma_prop = mvnrnd(gamma_mean_prop,gamma_cov_prop);
vari_new = exp(gamma_prop(1)+gamma_prop(2)*(Sa));
accept_ratio = sum(log(normpdf(log(Dr),(beta_prop(ii,:)*X')',sqrt(vari_new)))-(log(normpdf(log(Dr),(beta_prop(ii,:)*X')',sqrt(vari_old)))))+log(mvnpdf(gamma_prop',gamma_mean_pri,gamma_cov_pri))-log(mvnpdf(gamma_old',gamma_mean_pri,gamma_cov_pri));
decision = log(rand);
if decision<accept_ratio
gamma_new(ii,:) = gamma_prop;
vari_old = vari_new;
gamma_old = gamma_prop;
count = count+1;
%gamma_cov_pri = gamma_cov_prop;
%gamma_mean_pri = gamma_mean_prop;
else
gamma_new(ii,:) = gamma_old;    
end
end
parameters_bayes(1) = mean(beta_prop(1000:max_iter,1));
parameters_bayes(2) = mean(beta_prop(1000:max_iter,2));
parameters_bayes(3) = mean(gamma_new(1000:max_iter,1));
parameters_bayes(4) = mean(gamma_new(1000:max_iter,2));
%parameters_bayes(5) = mean(gamma_new(:,3));
cov_matrix_bayes(1:2,1:2) = cov(beta_prop(1000:max_iter,1),beta_prop(1000:max_iter,2));
cov_matrix_bayes(3:4,3:4) = cov(gamma_new(1000:max_iter,1),gamma_new(1000:max_iter,2));
%K = cov(gamma_new(:,2),gamma_new(:,3));
% cov_matrix_bayes(4,5) = K(1,2);
% cov_matrix_bayes(5,4) = K(1,2);
% cov_matrix_bayes(5,5) = K(2,2);
parameters_bayes'
cov_matrix_bayes
% hold on
% scatter(Sa,sqrt(K),'*')
sa_range = 0.05:0.05:2;
exact = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA4\Exact_stds_2g.txt');
reg_para = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA4\reg.txt');
plot(sa_range,sqrt((exp(parameters(3)+parameters(4)*(sa_range)))),'red','linewidth',1.5)
hold on
plot(sa_range,sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(sa_range))),sa_range,sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2))),'green','linewidth',1.5)%
hold on
% ex = interp1(exact(:,1),exact(:,2),sa_range);
% scatter(exact(:,1),exact(:,2))
% hold on
scatter(exact(:,1),exact(:,2))
%rBayes = rsquare(exact(:,2),sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(exact(:,1))+parameters_bayes(5)*(exact(:,1)).^2)));
SSE_Bayes = sum((exact(:,2)-sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(exact(:,1))))).^2);
%rFreq = rsquare(exact(:,2),sqrt((exp(parameters(3)+parameters(4)*(exact(:,1))+parameters(5)*(exact(:,1)).^2))));
SSE_Freq = sum((exact(:,2)-sqrt((exp(parameters(3)+parameters(4)*(exact(:,1)))))).^2);
%rFit = rsquare(exact(:,2),sqrt((exp(-3.83+2.91*(exact(:,1))-1.02*(exact(:,1)).^2))));
SSE_Fit = sum((exact(:,2)'-sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2)))).^2);
% rBayes
% rFreq
% rFit
std_freq = sqrt((exp(parameters(3)+parameters(4)*(sa_range))));
std_bayes = sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(sa_range)));
std_exact = sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2)));

SSE_Bayes
SSE_Freq
SSE_Fit
reg_para
grid on
xlabel('Sa (g)')
ylabel('Standard deviation in roof drift')
legend(strcat('Frequentist SSE = ',num2str(SSE_Freq)),strcat('Bayesian SSE = ',num2str(SSE_Bayes)),strcat('Fit to actual data SSE = ',num2str(SSE_Fit)))
%% Fisher's scoring algorithm
% Written by Somayajulu on March 25, 2016
% Modified by Somayajulu on April 1, 2016 
% CHANGE THE PATH IF NECESSARY
% Dr - Dependent, Sa - Independent
% NOTE: THE FIRST PART OF THE CODE IS AITKIN's ORIGINAL FISHER SCORING ALGORITHM
% NOTE: THE SECOND PART OF THE CODE IS GIBBS-METROPOLIS ALGORITHM WITH
% CEPEDA-GAMERMAN CORRECTION FOR DEPENDENT VARIABLE TO ENSURE GOOD
% ACCEPTANCE RATES FOR PROPOSALS (this is Bayesian)
% QUADRATIC VARIANCE FUNCTION 
clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\Roof Drift\';
Sa = importdata(strcat(file_path,'x.txt'));
%Sa = Sa(1:41);
%Sa = Sa(find(Sa<1));
Dr = importdata(strcat(file_path,'y.txt'));
%Dr = Dr(1:41);
%Dr = Dr(find(Sa<1));
X = [ones(max(size(Sa)),1) log(Sa) Sa];
Z = [ones(max(size(Sa)),1) (Sa) (Sa).^2];
beta_0 = regress(log(Dr),X);
res = log(Dr)-(beta_0(1)+beta_0(2)*log(Sa)+beta_0(3)*(Sa));
OLS_res = res;
vari = std(res)^2*ones(max(size(Sa)),1);
OLS_std = sqrt(vari);
gamma_0 = [log(std(res)^2);0;0];
parameters = vertcat(beta_0,gamma_0);
W = diag(1./vari);
A_0 = blkdiag(X'*W*X,0.5*Z'*Z);
temp = X'*(res./vari);
temp1 = 0.5*Z'*(res.^2./vari-1);
d_0 = vertcat(temp,temp1);
count = 1;
delta = 3;
while norm(delta)>=0.01
delta = inv(A_0)*d_0;
parameters = parameters+delta;
res = log(Dr)-(parameters(1)+parameters(2)*log(Sa)+parameters(3)*(Sa));
vari = (exp(parameters(3)+parameters(4)*(Sa)+parameters(5)*(Sa).^2));
temp = X'*(res./vari);
temp1 = 0.5*Z'*(res.^2./vari-1);
d_0 = vertcat(temp,temp1);
W = diag(1./vari);
A_0 = blkdiag(X'*W*X,0.5*Z'*Z);
count = count+1;
if count>1000
break
end
end
parameters_Freq = parameters;
parameters_Freq
%%
clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\RD\';
Sa = importdata(strcat(file_path,'x_pga.txt'));
% display('Independent is Sa3')
% fp = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\';
% eps = importdata(strcat(fp,'eps_sa1.txt'));
% display('eps is Sa3')
%Sa = Sa(1:41);
Dr = importdata(strcat(file_path,'y.txt'));
M = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Mw.txt'));
R = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Rjb.txt'));
eps = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\eps_pga.txt');
R(35) = 0.1;
R(36) = 0.1;
% M = importdata(strcat(fp,'Mw.txt'));
% R = importdata(strcat(fp,'Rjb.txt'));
X = [ones(max(size(Sa)),1) log(Sa) eps];% (R) eps];
Z = [ones(max(size(Sa)),1) (Sa) (Sa).^2];
beta_0 = regress(log(Dr),X);
res = log(Dr)-(beta_0(1)+beta_0(2)*log(Sa)+beta_0(3)*eps);%+beta_0(4)*(R));%+beta_0(4)*(R)+beta_0(5)*eps);
OLS_res = res;
vari = std(res);
%cov = regstats(log(Dr),X,'linear',{'covb'});
%cov.covb
str_eps = regstats(log(Dr),X,'linear',{'tstat'})
% str = regstats(res,M,'linear',{'tstat'});
% str1 = regstats(res,(R),'linear',{'tstat'});
% str2 = regstats(res,eps,'linear',{'tstat'});
% M1 = distcorr(M,res)
% R1 = distcorr(R,res)
% eps1 = distcorr(eps,res)
% final = [M1 R1 eps1];
colldiag(X)
cov_matrix = vari^2*inv(X'*X)