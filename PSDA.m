%% Roof Drift
clr
Sa = 0.001:0.001:2;
std_homo = ones(1,2000)*0.2825;
std_dev_exact = sqrt(exp(-2.938*Sa.^2+4.763*Sa-4.153));
std_bayes = sqrt(exp(-2.3816*Sa.^2+3.2921*Sa-3.3982));%-2.5457*Sa.^2+4.7817*Sa-4.424
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.8249*log(Sa)-3.447))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.001;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7659*log(Sa)-3.5008))./std_bayes,0,1);% 0.5977*log(Sa)-3.5394
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.001;

prob_homo = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7538*log(Sa)-3.4724))./std_homo,0,1);% 0.7568*log(Sa)-3.4777
AFE_homo(ii) = sum(prob_homo.*dsa1)*0.001;
end
% AFE_homo = -log(1-AFE_homo)/50;
% AFE_Bayes = -log(1-AFE_Bayes)/50;
% AFE_exact = -log(1-AFE_exact)/50;
AFE_exact = AFE_exact';
AFE_homo = AFE_homo';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_homo,'red','linewidth',1.5)
xlabel('Roof Drift')
%ylabel('Probability of exceedance in 50 years')
%plot(0.01:0.01:2,prob_Bayes,['red','--'],0.01:0.01:2,prob_exact,'green',0.01:0.01:2,prob_homo,['blue',':'],'linewidth',4)
legend('Exact','Bayes','Homo')
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 25)
  xlabel('Roof Drift')
  %ylabel('Pr.(roof drift > 0.04|IM)')  
  grid on
%% temp
clr
Sa = 0.01:0.01:1.5;
const_sd = 0.2681*ones(150,1);
bayes_sd = sqrt(exp(-2.0438*Sa.^2+3.1525*Sa-3.6207));%-2.0851*Sa.^2+3.944*Sa-4.3
std_dev_exact = sqrt(exp(-3.138*Sa.^2+4.861*Sa-3.846));%
std_dev_48 = sqrt(exp(-4.2427*Sa.^2+6.0394*Sa-4.4204));% -4.712*Sa.^2+7.2823*Sa-4.5567
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.839*log(Sa)-3.272))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.76*log(Sa)-3.4691))./bayes_sd,0,1);% 0.7726*log(Sa)-3.3639
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7605*log(Sa)-3.4687))./std_dev_48,0,1);% 0.8196*log(Sa)-3.3277
AFE_var(ii) = sum(prob_var.*dsa1)*0.01;

prob_const = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7533*log(Sa)-3.4570))./const_sd',0,1);
AFE_const(ii) = sum(prob_const.*dsa1)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
% AFE_Bayes = -log(1-AFE_Bayes)/50;
% AFE_exact = -log(1-AFE_exact)/50;
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red',dr_range,AFE_const,'m','linewidth',1.5)
xlabel('IDR 1')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Bayesian','Frequentist','Homoskedastic')
K = [dr_range;AFE_exact';AFE_Bayes';AFE_var';AFE_const];
%% temp1
clr
Sa = 0.01:0.01:2;
const_sd = 0.32*ones(200,1);
bayes_sd = sqrt(exp(-0.558*Sa-2.28));%-2.0172*Sa.^2+3.1211*Sa-3.611
std_dev_exact = sqrt(exp(-0.379*Sa.^2+0.187*Sa-2.251));
std_dev_48 = sqrt(exp(-0.5959*Sa-2.305));% -4.712*Sa.^2+7.2823*Sa-4.5567
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 1:1500;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.8257*log(Sa)+6.01))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7654*log(Sa)+6.0727))./bayes_sd,0,1);% 0.7726*log(Sa)-3.3639
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7663*log(Sa)+6.0733))./std_dev_48,0,1);% 0.8196*log(Sa)-3.3277
AFE_var(ii) = sum(prob_var.*dsa1)*0.01;

prob_const = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.779*log(Sa)+6.0805))./const_sd',0,1);
AFE_const(ii) = sum(prob_const.*dsa1)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
% AFE_Bayes = -log(1-AFE_Bayes)/50;
% AFE_exact = -log(1-AFE_exact)/50;
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red',dr_range,AFE_const,'m','linewidth',1.5)
xlabel('PFA1 (in/sec^2)')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Bayesian','Frequentist','Homoskedastic')
K = [dr_range;AFE_exact';AFE_Bayes';AFE_var';AFE_const];
%% IDR 1
clr
Sa = 0.01:0.01:1.5;
bayes_sd = ones(1,150)*0.3137;%0.3704;
std_dev_exact = sqrt(exp(-3.138*Sa.^2+4.861*Sa-3.846));
std_dev_48 = sqrt(exp(-4.7121*Sa.^2+7.2823*Sa-4.5567));% -4.712*Sa.^2+7.2823*Sa-4.5567
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.9289*log(Sa)-3.31))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.8656*log(Sa)-3.31))./bayes_sd,0,1);% 0.7726*log(Sa)-3.3639
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.8182*log(Sa)-3.3879))./std_dev_48,0,1);% 0.8196*log(Sa)-3.3277
AFE_var(ii) = sum(prob_var.*dsa1)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red','linewidth',1.5)
xlabel('IDR 1')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Homoskedastic','Heteroskedastic')
%% IDR 2
clr
Sa = 0.01:0.01:1.8;
bayes_sd = ones(1,180)*0.4126;
std_dev_exact = sqrt(exp(-0.1793*Sa.^2-0.01961*Sa-1.546));
std_dev_48 = sqrt(exp(-1.5876*Sa.^2+1.4932*Sa-2.033));%-3.5934*Sa.^2+6.525*Sa-4.988
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.767*log(Sa)-3.617))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.4917*log(Sa)-3.6770))./bayes_sd,0,1);% 0.7388*log(Sa)-3.3371
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.4887*log(Sa)-3.6817))./std_dev_48,0,1);% 0.7962*log(Sa)-3.2706
AFE_var(ii) = sum(prob_var.*dsa1)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red','linewidth',1.5)
xlabel('IDR2')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Homoskedastic','Heteroskedastic')
%% IDR 3
clr
Sa = 0.01:0.01:1.8;
bayes_sd = ones(1,180)*0.3442;
std_dev_exact = sqrt(exp(0.047*Sa.^2-0.5257*Sa-1.72));
std_dev_48 = sqrt(exp(-3.5763*Sa.^2+2.8556*Sa-2.4618));%-3.5935*Sa.^2+6.5254*Sa-4.9883
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7419*log(Sa)-3.598))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.5224*log(Sa)-3.6335))./bayes_sd,0,1);% 0.7388*log(Sa)-3.3371
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.4887*log(Sa)-3.6766))./std_dev_48,0,1);% 0.7962*log(Sa)-3.2706
AFE_var(ii) = sum(prob_var.*dsa1)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red','linewidth',1.5)
xlabel('IDR2')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Homoskedastic','Heteroskedastic')
%% IDR 4
clr
Sa = 0.01:0.01:1.5;
bayes_sd = ones(1,150)*0.3023;
std_dev_exact = sqrt(exp(-1.138*Sa.^2+1.478*Sa-2.63));
%std_dev_48 = sqrt(exp(-3.5925*Sa.^2+5.8174*Sa-4.4154));%-3.5935*Sa.^2+6.5254*Sa-4.9883
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.6553*log(Sa)-3.567))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.42183*log(Sa)-3.703928))./bayes_sd,0,1);% 0.7388*log(Sa)-3.3371
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.01;

% prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.6674*log(Sa)-3.3340))./std_dev_48,0,1);% 0.7962*log(Sa)-3.2706
% AFE_var(ii) = sum(prob_var.*dsa1)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
AFE_exact = AFE_exact';
%AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green','linewidth',1.5);%,dr_range,AFE_var,'red'
xlabel('IDR2')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Homoskedastic','Heteroskedastic')
%% Joint Rotation
clr
Sa = 0.01:0.01:1.5;
bayes_sd = ones(1,150)*0.4624;
% 80 points mean:0.8361*log(Sa)-3.4143   sd:0.3040
% 40 points mean:0.8769*log(Sa)-3.39853   sd:0.31932
std_dev_exact = sqrt(exp(-0.2134*Sa.^2-0.01323*Sa-1.33));
std_dev_48 = sqrt(exp(-1.9417*Sa.^2+1.2584*Sa-1.5831));
% 80 points mean:0.8941*log(Sa)-3.3629   sd:-3.222*Sa.^2+5.1223*Sa-3.9467
% 40 points mean:0.9822*log(Sa)-3.2705   sd:-4.9077*Sa.^2+8.4938*Sa-5.2277
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 0.001:0.001:0.08;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.8579*log(Sa)-3.748))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.5912*log(Sa)-3.787))./bayes_sd,0,1);
AFE_Bayes(ii) = sum(prob_Bayes.*dsa1)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.5834*log(Sa)-3.7963))./std_dev_48,0,1);
AFE_var(ii) = sum(prob_var.*dsa1)*0.01;
end
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red','linewidth',1.5)
xlabel('Joint Rotation')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Homoskedastic','Heteroskedastic')
%% PFA 1
clr
Sa = 0.01:0.01:1.8;
bayes_sd = ones(1,180)*0.2655;
std_dev_exact = sqrt(exp(-0.1636*Sa.^2+0.01377*Sa-2.3));
std_dev_48 = sqrt(exp(0.5971*Sa.^2-0.7720*Sa-2.4765));%0.59741*Sa.^2-0.77202*Sa-2.47649
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa = abs((Differentiation((0.01),AFE_sa1)));
dr_range = 1:1:1500;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.8257*log(Sa)+6.01))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7635*log(Sa)+6.0519))./bayes_sd,0,1);% 0.74889*log(Sa)+6.038
AFE_Bayes(ii) = sum(prob_Bayes.*dsa)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.7489*log(Sa)+6.0383))./std_dev_48,0,1);% 0.7634*log(Sa)+6.0518
AFE_var(ii) = sum(prob_var.*dsa)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red','linewidth',1.5)
xlabel('PFA1 (cm/sec^2)')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Homoskedastic','Heteroskedastic')
%% PFA 2
clr
Sa = 0.01:0.01:2;
bayes_sd = ones(1,200)*0.179911;
std_dev_exact = sqrt(exp(0.2962*Sa.^2-0.9652*Sa-1.047));
std_dev_48 = sqrt(exp(-0.72812*Sa.^2-0.10896*Sa-3.139703));%0.59741*Sa.^2-0.77202*Sa-2.47649
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
pga_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_pga = exp(interp1(log(pga_temp(:,1)),log(pga_temp(:,2)),log(Sa),'spline'));
dpga = abs((Differentiation((0.01),AFE_pga)));
dr_range = 1:1:1500;
for ii = 1:max(size(dr_range))
prob_exact = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.737*log(Sa)+6.077))./std_dev_exact,0,1);
AFE_exact(ii) = sum(prob_exact.*dsa1)*0.01;

prob_Bayes = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.794057*log(Sa)+6.2666))./bayes_sd,0,1);% 0.2095*log(Sa)+5.5675
AFE_Bayes(ii) = sum(prob_Bayes.*dpga)*0.01;

prob_var = 1-normcdf((log(dr_range(ii))*ones(max(size(Sa)),1)'-(0.806756*log(Sa)+6.27716))./std_dev_48,0,1);% 0.7962*log(Sa)-3.2706
AFE_var(ii) = sum(prob_var.*dpga)*0.01;
end
% AFE_const = -log(1-AFE_const)/50;
% AFE_var = -log(1-AFE_var)/50;
AFE_exact = AFE_exact';
AFE_var = AFE_var';
AFE_Bayes = AFE_Bayes';
semilogy(dr_range,AFE_exact,'blue',dr_range,AFE_Bayes,'green',dr_range,AFE_var,'red','linewidth',1.5)
xlabel('PFA1 (cm/sec^2)')
ylabel('Probability of exceedance in 50 years')
legend('Exact','Homoskedastic','Heteroskedastic')