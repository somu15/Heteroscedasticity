clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\PFA4\';
name = 'pgv';
Sa = importdata(strcat(file_path,'x_',name,'.txt'));
eps = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\eps_',name,'.txt'));

Dr = importdata(strcat(file_path,'y.txt'));
M = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Mw.txt'));
R = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Rjb.txt'));
Vs = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Vs30.txt'));

X = [ones(max(size(Sa)),1) log(Sa)];
beta = regress(log(Dr),X);
res = log(Dr)-(beta(1)+beta(2)*log(Sa));
pear_M = corr(M,res)
pear_R = corr(R,res)
pear_eps = corr(eps,res)
dist_M = distcorr(M,res)
dist_R = distcorr(R,res)
dist_eps = distcorr(eps,res)
dist = [dist_M dist_R dist_eps];
pear = [pear_M pear_R pear_eps];