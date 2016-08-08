%% Function for the Fisher scoring and Bayesian algorithms
%
% Capturing heteroscedasticity from Dhulipala and Flint, 2016.
% ""
%
% Created by Somayajulu Dhulipala on March 25, 2016
% Modified by Somayajulu Dhulipala on April 1, 2016 
% Modified by Somayajulu Dhulipala on July 29, 2016

% Main references
% [1] "Modelling Variance Heterogeneity in Normal Regression using GLIM" 
%     by Murray Aitkin, 1987, in Journal of the Royal Stat. Soc., series C
%     (Applied Statistics).
% [2] "Bayesian Modeling of Variance Heterogeneity in Normal Regression
%     Models" by Edilberto Cepeda and Dani Gamerman, 2001, in Brazilian 
%     Journal of Probability and Statistics.




function [parameters_Freq] = Fisher_Bayes(IM,Response)
%% Create variables for OLS regression
X = [ones(max(size(IM)),1) log(IM)];
Z = [ones(max(size(IM)),1) (IM) (IM).^2];
%% Perform OLS and compute the constant variance vector
beta_0 = regress(log(Response),X);
res = log(Response)-(beta_0(1)+beta_0(2)*log(IM));
vari = std(res)^2*ones(max(size(IM)),1);
%% Create more variables for the Fisher scoring algorithm
gamma_0 = [log(std(res)^2);0;0];
parameters = vertcat(beta_0,gamma_0);
%% Perform the first iteration in the Fisher scoring algorithm
W = diag(1./vari);
A_0 = blkdiag(X'*W*X,0.5*(Z'*Z));
temp = X'*(res./vari);
temp1 = 0.5*Z'*(res.^2./vari-1);
d_0 = vertcat(temp,temp1);
%% Specify the stopping criteria and associated tracking variables
count = 1;
stop_crit_count = 1000;
delta = 3; % This the norm of the error vector which needs to be minimized. 
% Initial guess can be arbitrarily greater than the stopping norm.
stop_crit_delta = 0.01;
%% Fisher scoring algorithm
while norm(delta) >= stop_crit_delta
    delta = inv(A_0)*d_0;
    parameters = parameters+delta;
    res = log(Response)-(parameters(1)+parameters(2)*log(IM));
    vari = (exp(parameters(3)+parameters(4)*(IM)+parameters(5)*(IM).^2));
    temp = X'*(res./vari);
    temp1 = 0.5*Z'*(res.^2./vari-1);
    d_0 = vertcat(temp,temp1);
    W = diag(1./vari);
    A_0 = blkdiag(X'*W*X,0.5*(Z'*Z));
    count = count+1;
    if count > stop_crit_count
        break
    end
end
parameters_Freq = parameters;
%% 
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
beta_mean_post = beta_cov_post*(inv(beta_cov_pri)*beta_mean_pri+X'*diag(vari_old)*log(Response));
beta_prop(ii,:) = mvnrnd(beta_mean_post,beta_cov_post);
Trans = log(vari_old)+(log(Response)-(beta_prop(ii,:)*X')').^2./vari_old-1;
gamma_cov_prop = inv(inv(gamma_cov_pri)+0.5*Z'*Z);
gamma_mean_prop = gamma_cov_prop*(inv(gamma_cov_pri)*gamma_mean_pri+0.5*Z'*Trans);
gamma_prop = mvnrnd(gamma_mean_prop,gamma_cov_prop);
vari_new = exp(gamma_prop(1)+gamma_prop(2)*(IM)+gamma_prop(3)*(IM).^2);
accept_ratio = sum(log(normpdf(log(Response),(beta_prop(ii,:)*X')',sqrt(vari_new)))-(log(normpdf(log(Response),(beta_prop(ii,:)*X')',sqrt(vari_old)))))+log(mvnpdf(gamma_prop',gamma_mean_pri,gamma_cov_pri))-log(mvnpdf(gamma_old',gamma_mean_pri,gamma_cov_pri));
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
subplot(1,2,2)
hist(IM)
xlabel('Sa(T1)')
ylabel('Frequency')
subplot(1,2,1)
sa_range = 0.05:0.05:2;
exact = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA1\Exact_stds_2g.txt');
reg_para = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\PFA1\reg.txt');
plot(sa_range,sqrt((exp(parameters(3)+parameters(4)*(sa_range)+parameters(5)*(sa_range).^2))),'red','linewidth',1.5)
hold on
plot(sa_range,sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(sa_range)+parameters_bayes(5)*(sa_range).^2)),sa_range,sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2))),'green','linewidth',1.5)%
hold on
scatter(exact(:,1),exact(:,2))
SSE_Bayes = sum((exact(:,2)-sqrt(exp(parameters_bayes(3)+parameters_bayes(4)*(exact(:,1))+parameters_bayes(5)*(exact(:,1)).^2))).^2);
SSE_Freq = sum((exact(:,2)-sqrt((exp(parameters(3)+parameters(4)*(exact(:,1))+parameters(5)*(exact(:,1)).^2)))).^2);
SSE_Fit = sum((exact(:,2)'-sqrt((exp(reg_para(1)+reg_para(2)*(sa_range)+reg_para(3)*(sa_range).^2)))).^2);
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
end