%% Fisher's scoring algorithm
% Written by Somayajulu on March 25, 2016
% CHANGE THE PATH IF NECESSARY
% Dr - Dependent, Sa - Independent
% NOTE: THE FIRST PART OF THE CODE IS AITKIN's ORIGINAL FISHER SCORING ALGORITHM
% NOTE: THE SECOND PART OF THE CODE IS GIBBS-METROPOLIS ALGORITHM WITH
% CEPEDA-GAMERMAN CORRECTION FOR DEPENDENT VARIABLE TO ENSURE GOOD
% ACCEPTANCE RATES FOR PROPOSALS
% QUADRATIC VARIANCE FUNCTION 
%clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\';
Sa = importdata(strcat(file_path,'Sa_more_points.txt'));
%Sa = Sa(40:80);
Dr = importdata(strcat(file_path,'Drift_more_points.txt'));
%Dr = Dr(40:80);
X = [ones(max(size(Sa)),1) log(Sa)];
Z = [ones(max(size(Sa)),1) (Sa) (Sa).^2];
beta_0 = regress(log(Dr),X);
res = log(Dr)-(beta_0(1)+beta_0(2)*log(Sa));
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
res = log(Dr)-(parameters(1)+parameters(2)*log(Sa));
vari = (exp(parameters(3)+parameters(4)*(Sa)+parameters(5)*(Sa).^2));
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
gamma_mean_pri = 0.0001*parameters(3:max(size(parameters)));
cov_matrix = inv(A_0);
beta_cov_pri = cov_matrix(1:2,1:2);
gamma_cov_pri = cov_matrix(3:max(size(parameters)),3:max(size(parameters)));
gamma_cov_pri = 10^10*gamma_cov_pri;
gamma_old = mvnrnd(gamma_mean_pri,gamma_cov_pri);
vari_old = vari;
count = 0;
max_iter = 1000;
for ii = 1:max_iter
beta_cov_post = inv(inv(beta_cov_pri)+X'*diag(vari_old)*X);
beta_mean_post = beta_cov_post*(inv(beta_cov_pri)*beta_mean_pri+X'*diag(vari_old)*log(Dr));
beta_new(ii,:) = mvnrnd(beta_mean_post,beta_cov_post);
Trans = log(vari_old)+(log(Dr)-(beta_new(ii,:)*X')').^2./vari_old-1;
gamma_cov_prop = inv(inv(gamma_cov_pri)+0.5*Z'*Z);
gamma_mean_prop = gamma_cov_prop*(inv(gamma_cov_pri)*gamma_mean_pri+0.5*Z'*Trans);
gamma_prop = mvnrnd(gamma_mean_prop,gamma_cov_prop);
vari_new = exp(gamma_prop(1)+gamma_prop(2)*(Sa)+gamma_prop(3)*(Sa).^2);
accept_ratio = sum(log(normpdf(log(Dr),(beta_new(ii,:)*X')',sqrt(vari_new)))-(log(normpdf(log(Dr),(beta_new(ii,:)*X')',sqrt(vari_old)))))+log(mvnpdf(gamma_prop',gamma_mean_pri,gamma_cov_pri))-log(mvnpdf(gamma_old',gamma_mean_pri,gamma_cov_pri));
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
parameters_bayes(1) = mean(beta_new(:,1));
parameters_bayes(2) = mean(beta_new(:,2));
parameters_bayes(3) = mean(gamma_new(100:max_iter,1));
parameters_bayes(4) = mean(gamma_new(100:max_iter,2));
parameters_bayes(5) = mean(gamma_new(100:max_iter,3));
cov_matrix_bayes(1:2,1:2) = cov(beta_new(100:max_iter,1),beta_new(100:max_iter,2));
cov_matrix_bayes(3:4,3:4) = cov(gamma_new(100:max_iter,1),gamma_new(100:max_iter,2));
K = cov(gamma_new(100:max_iter,2),gamma_new(100:max_iter,3));
cov_matrix_bayes(4,5) = K(1,2);
cov_matrix_bayes(5,4) = K(1,2);
cov_matrix_bayes(5,5) = K(2,2);
parameters_bayes'
cov_matrix_bayes
K = exp(parameters_bayes(3)+parameters_bayes(4)*(Sa)+parameters_bayes(5)*(Sa).^2);
figure()
subplot(3,2,1)
plot(10:max_iter,beta_new(10:max_iter,1))
grid on
xlabel('Iteration')
ylabel('Beta 1')
subplot(3,2,2)
plot(10:max_iter,beta_new(10:max_iter,2))
grid on
xlabel('Iteration')
ylabel('Beta 2')
subplot(3,2,3)
plot(10:max_iter,gamma_new(10:max_iter,1))
grid on
xlabel('Iteration')
ylabel('Gamma 1')
subplot(3,2,4)
plot(10:max_iter,gamma_new(10:max_iter,2))
grid on
xlabel('Iteration')
ylabel('Gamma 2')
subplot(3,2,5)
plot(10:max_iter,gamma_new(10:max_iter,3))
grid on
xlabel('Iteration')
ylabel('Gamma 3')
subplot(3,2,6)
scatter((Sa),sqrt(vari))
hold on
scatter((Sa),sqrt(K),'*')
grid on
xlabel('Sa (g)')
ylabel('Standard deviation in roof drift')
legend('Fisher scoring algorithm','Fisher-Gibbs-Metropolis algorithm')
%% Fisher's scoring algorithm
% Written by Somayajulu on March 25, 2016
% CHANGE THE PATH IF NECESSARY
% Dr - Dependent, Sa - Independent
% NOTE: THE FIRST PART OF THE CODE IS AITKIN's ORIGINAL FISHER SCORING ALGORITHM
% NOTE: THE SECOND PART OF THE CODE IS GIBBS-METROPOLIS ALGORITHM WITH
% CEPEDA-GAMERMAN CORRECTION FOR DEPENDENT VARIABLE TO ENSURE GOOD
% ACCEPTANCE RATES FOR PROPOSALS
% QUARTIC VARIANCE FUNCTION 
clear all
clc
file_path  = 'C:\Users\SOMAYAJULU\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\';
Sa = importdata(strcat(file_path,'Sa_more_points.txt'));
%Sa = Sa(40:80);
Dr = importdata(strcat(file_path,'Drift_more_points.txt'));
%Dr = Dr(40:80);
X = [ones(max(size(Sa)),1) log(Sa)];
Z = [ones(max(size(Sa)),1) (Sa) (Sa).^2 Sa.^3 Sa.^4];
beta_0 = regress(log(Dr),X);
res = log(Dr)-(beta_0(1)+beta_0(2)*log(Sa));
OLS_res = res;
vari = std(res)^2*ones(max(size(Sa)),1);
OLS_std = sqrt(vari);
gamma_0 = [log(std(res)^2);0;0;0;0];
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
vari = (exp(parameters(3)+parameters(4)*(Sa)+parameters(5)*(Sa).^2+parameters(6)*(Sa).^3+parameters(7)*(Sa).^4));
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
cov_matrix = inv(A_0);
cov_matrix