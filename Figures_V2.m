%% Figure 1
Sa = 0.01:0.01:2;
prob_ex_Sa = 1-normcdf((log(300)*ones(max(size(Sa)),1)'-(0.6251*log(Sa)+6.0115))./(0.3*ones(1,max(size(Sa)))),0,1);
prob_ex1 = 1-normcdf((log(300)*ones(max(size(Sa)),1)'-(0.6251*log(Sa)+6.0115))./(0.00001*ones(1,max(size(Sa)))),0,1);
subplot(1,2,1)
c1 = [0.3010 0.745 0.933];
c2 = [0.635 0.078 0.184];
plot(Sa,prob_ex1,'--','color',c1,'linewidth',3)
hold on
plot(Sa,prob_ex_Sa,'color',c2,'linewidth',3)
hold on
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')  
set(gca, 'FontSize', 20)
  xlabel('IM')
  ylabel('Pr.(\theta_{max} > x|IM)')
  title('(a)')
plot(0.61,0.5,['.','g'],'MarkerSize',40)
legend('Zero standard deviation','Non-zero standard deviation','50% exceedance probability')
prob_ex_Sa = 1-normcdf((log(300)*ones(max(size(Sa)),1)'-(0.6251*log(Sa)+6.0115))./(0.3*ones(1,max(size(Sa)))),0,1);
prob_ex1 = 1-normcdf((log(300)*ones(max(size(Sa)),1)'-(0.4*log(Sa)+6.0115))./(0.2*ones(1,max(size(Sa)))),0,1);
subplot(1,2,2)
c1 = [0.3010 0.745 0.933];
c2 = [0.635 0.078 0.184];
plot(Sa,prob_ex1,'--','color',c1,'linewidth',3)
hold on
plot(Sa,prob_ex_Sa,'color',c2,'linewidth',3)
hold on
set(gca,'YTick',[])
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 20)
  xlabel('IM')
  title('(b)')
  legend('Without conditional independence','With conditional independence')
  %% Figure 2
  Sa = 0.01:0.01:2;
  dIM = 0.01;
  sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
prob_ex_Sa = 1-normcdf((log(300)*ones(max(size(Sa)),1)'-(0.6251*log(Sa)+6.0115))./(0.0000001*ones(1,max(size(Sa)))),0,1);
dsa1 = abs((Differentiation((dIM),AFE_sa1)));
subplot(1,2,1)
c1 = [0.929 0.694 0.125];
c2 = [0.85 0.325 0.098];
c1 = [0.3010 0.745 0.933];
c2 = [0.635 0.078 0.184];
loglog(Sa,AFE_sa1,'color',c1,'linewidth',3)
set(gca,'YTick',[])
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 25)
  xlabel('IM*')
  ylabel('\lambda(IM*)')
  title('(a)')
  subplot(1,2,2)
  bayes_Sa = prob_ex_Sa.*dsa1./sum(prob_ex_Sa.*dsa1*dIM);
bayes_Sa = bayes_Sa/trapz(Sa,bayes_Sa);
% [AX,H1,H2] = plotyy(x,y1,x,y2,'plot');
% set(AX,{'ycolor'},{'r';'b'})  % Left color red, right color blue...
yyaxis left
[haxes] = plot(Sa,prob_ex_Sa,'color',c1,'linewidth',3);
set(gca, 'FontSize', 20)
set(gca,'YTick',[])
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')
 xlabel('IM*');
t = ylabel('Pr.(\theta_{max} > x |IM*)');
  t.Color = 'k';
  yyaxis right
  plot(Sa,bayes_Sa,'--','color',c2,'linewidth',3)
  set(gca, 'FontSize', 20)
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')
 t =  ylabel('w(IM*)');
 t.Color = 'k';
  title('(b)')
  hold on
hold on
  plot(0.61,0.1,['.','g'],'MarkerSize',40)
  legend('Demand fragility','Weight function','IMc*')
  %% Figure 3
  Sa = 0.01:0.01:2;
  dIM = 0.01;
  sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
prob_ex_Sa = 1-normcdf((log(300)*ones(max(size(Sa)),1)'-(0.6251*log(Sa)+6.0115))./(0.25*ones(1,max(size(Sa)))),0,1);
dsa1 = abs((Differentiation((dIM),AFE_sa1)));
 bayes_Sa = prob_ex_Sa.*dsa1./sum(prob_ex_Sa.*dsa1*dIM);
bayes_Sa = bayes_Sa/trapz(Sa,bayes_Sa);
     subplot(1,2,1)
     c1 = [0.85 0.3250 0.098];
   loglog(Sa,AFE_sa1,'color',c1,'linewidth',4)
   set(gca, 'FontName', 'Times')
   set(gca, 'FontSize', 20)
   xlabel('IM*')
   ylabel('F(IM*)')
   grid on
   subplot(1,2,2)
   yyaxis left
  
   plot(Sa,prob_ex_Sa,'color',c1,'linewidth',3)
   set(gca, 'FontName', 'Times')
   set(gca, 'FontSize', 20)
   xlabel('IM*')
   ylabel('Pr.(theta_{max} > x | IM*)')
   grid on
   yyaxis right
   plot(Sa,bayes_Sa,['red','--'],'linewidth',4)
   ylabel('w*')
   legend('CDF','weight function')
   hold on
   scatter(0.22,0,'g','linewidth',10)
   
  subplot(3,2,1)
  c1 = [0.3010 0.745 0.933];
c2 = [0.635 0.078 0.184];
c1 = [0.3010 0.745 0.933];
c2 = [0.635 0.078 0.184];
   loglog(Sa,AFE_sa1,'color',c1,'linewidth',3)
   set(gca, 'FontSize', 20)
set(gca, 'FontName', 'Times')
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  %xlabel('IM')
  ylabel('\lambda(IM)')
  title('(a)')
  file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\Roof Drift\';
Sa1 = importdata(strcat(file_path,'xtemp.txt'));
Dr = importdata(strcat(file_path,'ytemp.txt'));
  subplot(3,2,2)
  plot(log(Sa),(0.6702*log(Sa)-3.6251),'color',c1,'linewidth',4)
  set(gca, 'FontSize', 20)
  hold on
  scatter(log(Sa1),log(Dr),'filled','markerfacecolor',[0.494 0.184 0.556],'linewidth',2)
  set(gca, 'FontName', 'Times')
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  %xlabel('ln(IM)')
  ylabel('ln(\theta_{max})')
  title('(b)')
  subplot(3,2,3)
  plot(Sa,prob_ex_Sa,'color',c1,'linewidth',3)
  set(gca, 'FontSize', 20)
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')
  xlabel('IM')
  ylabel('Pr.(\theta_{max} > x|IM)')
  title('(c)')
  subplot(3,2,4)
  plot(Sa,bayes_Sa,'color',c1,'linewidth',3)
  set(gca, 'FontSize', 20)
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')
  xlabel('IM')
  ylabel('w(IM)')
  title('(d)')
  subplot(3,2,[5,6])
  plot(Sa,prob_ex_Sa,'color',c1,'linewidth',3)
  hold on
  temp = zeros(100,1);
    temp1 = ones(max(size(Sa))-100,1);
    psuedo_Sa = vertcat(temp,temp1);
  plot(Sa,psuedo_Sa,'--','color',c2,'linewidth',3)
  set(gca, 'FontSize', 20)
  set(gca,'YTick',[])
  set(gca,'XTick',[])
  set(gca, 'FontName', 'Times')
  xlabel('IM')
  ylabel('Pr.(\theta_{max} > x|IM)')
  title('(e)')
  legend('Realistic IM','Psuedo IM')
  %% Figure 4
  x = 1:6;
  y = [3.312 3.382 4.443 2.736 3.526 3.680;0.407 0.414 0.449 0.361 0.095 0.165;5.344 5.010 5.223 5.213 0.496 1.109;2.009 1.963 3.374 1.389 2.179 2.208;1.681 1.782 2.779 1.219 2.358 1.992];
  spider(y','A metric variation',[6 6 6 6 6 6]',{'RD 0.01','IDR1 0.01','IDR4 0.01','JR 0.01','PFA1 0.63g','PFA4 0.88g'},{'PGA','PGV','SaT1','SaT2','SaT3'})
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)

  %% Figure 5
  SaT1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\Roof Drift\x.txt');
  PGA = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\Roof Drift\x_alt.txt');
  subplot(1,2,1)
  hist(SaT1,25)
  set(get(gca,'child'),'FaceColor',[51/255 51/255 255/255],'EdgeColor','k');
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 25)
  xlabel('SaT1')
  ylabel('Frequency')
  title('(a) 83 data points')
  subplot(1,2,2)
  hist(PGA,25)
  set(get(gca,'child'),'FaceColor',[0.466 0.674 0.188],'EdgeColor','k');
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 25)
  xlabel('PGA')
  ylabel('Frequency')
  title('(b) 79 data points')
  %% Figure 6
  dr = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\A_conv_drift.txt');
  acc = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\A_conv_acc.txt');
  DR = 0.001:0.001:0.04;
  ACC = (1:10:772)*2.589e-3;
  c1 = [0.466 0.674 0.188];
  c2 = [51/255 51/255 255/255];
  for ii = 1:min(size(dr))/2
  dr_err(ii) = sum(abs(dr(2*ii-1,:)-dr(2*ii,:))/dr(2*ii-1,:)*100);
  end
  for ii = 1:min(size(acc))/2
  acc_err(ii) = sum(abs(acc(2*ii-1,:)-acc(2*ii,:))/acc(2*ii-1,:)*100);
  end
  subplot(3,4,1)
  plot(DR,dr(1,:),'color',c2,'linewidth',2.5)
  hold on
  plot(DR,dr(2,:),':','color',c2,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
legend('SaT1 exact','SaT1 approx')
grid on
title(strcat('RD (',num2str(dr_err(1)),'%)'))

  subplot(3,4,2)
  plot(DR,dr(3,:),'color',c1,'linewidth',2.5)
  hold on
  plot(DR,dr(4,:),':','color',c1,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 10])
set(gca,'XTick',[0 0.02 0.04])
set(gca,'yTick',[0 5 10])
grid on
legend('PGA exact','PGA approx')
title(strcat('RD (',num2str(dr_err(2)),'%)'))

subplot(3,4,5)
plot(DR,dr(5,:),'color',c2,'linewidth',2.5)
  hold on
  plot(DR,dr(6,:),':','color',c2,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 10])
ylabel('A')
set(gca,'XTick',[0 0.02 0.04])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('IDR1 (',num2str(dr_err(3)),'%)'))

subplot(3,4,6)
plot(DR,dr(7,:),'color',c1,'linewidth',2.5)
  hold on
  plot(DR,dr(8,:),':','color',c1,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 10])
set(gca,'XTick',[0 0.02 0.04])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('IDR1 (',num2str(dr_err(4)),'%)'))

subplot(3,4,7)
plot(DR,dr(9,:),'color',c2,'linewidth',2.5)
  hold on
  plot(DR,dr(10,:),':','color',c2,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 10])
set(gca,'XTick',[0 0.02 0.04])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('IDR4 (',num2str(dr_err(5)),'%)'))

subplot(3,4,8)
plot(DR,dr(11,:),'color',c1,'linewidth',2.5)
  hold on
  plot(DR,dr(12,:),':','color',c1,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 10])
set(gca,'XTick',[0 0.02 0.04])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('IDR4 (',num2str(dr_err(6)),'%)'))

subplot(3,4,3)
plot(DR,dr(13,:),'color',c2,'linewidth',2.5)
  hold on
  plot(DR,dr(14,:),':','color',c2,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 10])
set(gca,'XTick',[0 0.02 0.04])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('JR (',num2str(dr_err(7)),'%)'))

subplot(3,4,4)
plot(DR,dr(15,:),'color',c1,'linewidth',2.5)
  hold on
  plot(DR,dr(16,:),':','color',c1,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 10])
set(gca,'XTick',[0 0.02 0.04])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('JR (',num2str(dr_err(8)),'%)'))

subplot(3,4,9)
plot(ACC,acc(1,:),'color',c2,'linewidth',2.5)
  hold on
  plot(ACC,acc(2,:),':','color',c2,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 2 0 10])
set(gca,'XTick',[0 1 2])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('PFA1 (g) (',num2str(acc_err(1)),'%)'))

subplot(3,4,10)
plot(ACC,acc(3,:),'color',c1,'linewidth',2.5)
  hold on
  plot(ACC,acc(4,:),':','color',c1,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 2 0 10])
set(gca,'XTick',[0 1 2])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('PFA1 (g) (',num2str(acc_err(2)),'%)'))

subplot(3,4,11)
plot(ACC,acc(5,:),'color',c2,'linewidth',2.5)
  hold on
  plot(ACC,acc(6,:),':','color',c2,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 2 0 10])
set(gca,'XTick',[0 1 2])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('PFA4 (g) (',num2str(acc_err(3)),'%)'))

subplot(3,4,12)
plot(ACC,acc(7,:),'color',c1,'linewidth',2.5)
  hold on
  plot(ACC,acc(8,:),':','color',c1,'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 2 0 10])
set(gca,'XTick',[0 1 2])
set(gca,'yTick',[0 5 10])
grid on
title(strcat('PFA4 (g) (',num2str(acc_err(4)),'%)'))
%% Figure 7
dr = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\Final_result_Dr.txt');
  acc = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\Final_result_acc.txt');
  DR = 0.001:0.001:0.04;
  ACC = (1:10:772)*2.589e-3;
  col = [ 0.466 0.674 0.188
    204/255 0 102/255
    51/255 51/255 255/255
    153/255 153/255 255/255
    51/255 255/255 255/255
    51/255 255/255 255/255
    51/255 255/255 255/255];
  subplot(3,2,1)
  plot(DR,dr(2,:),'--','color',col(1,:),'linewidth',2.5)
  hold on
  plot(DR,dr(1,:),'color',col(3,:),'linewidth',2.5)
  hold on
  plot(DR,dr(3,:),':','color',col(4,:),'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 8])
set(gca,'XTick',[0 0.01 0.02 0.03 0.04])
set(gca,'yTick',[0 2 4 6 8])
grid on
xlabel('Drift')
title('RD')

subplot(3,2,3)
  plot(DR,dr(5,:),'--','color',col(1,:),'linewidth',2.5)
  hold on
  plot(DR,dr(4,:),'color',col(3,:),'linewidth',2.5)
  hold on
  plot(DR,dr(6,:),':','color',col(4,:),'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 8])
set(gca,'XTick',[0 0.01 0.02 0.03 0.04])
set(gca,'yTick',[0 2 4 6 8])
grid on
xlabel('Drift')
ylabel('A')
title('IDR1')

subplot(3,2,4)
  plot(DR,dr(8,:),'--','color',col(1,:),'linewidth',2.5)
  hold on
  plot(DR,dr(7,:),'color',col(3,:),'linewidth',2.5)
  hold on
  plot(DR,dr(9,:),':','color',col(4,:),'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 8])
set(gca,'XTick',[0 0.01 0.02 0.03 0.04])
set(gca,'yTick',[0 2 4 6 8])
grid on
xlabel('Drift')
title('IDR4')

subplot(3,2,2)
  plot(DR,dr(11,:),'--','color',col(1,:),'linewidth',2.5)
  hold on
  plot(DR,dr(10,:),'color',col(3,:),'linewidth',2.5)
  hold on
  plot(DR,dr(12,:),':','color',col(4,:),'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 0.04 0 8])
set(gca,'XTick',[0 0.01 0.02 0.03 0.04])
set(gca,'yTick',[0 2 4 6 8])
grid on
xlabel('Drift')
title('JR')

subplot(3,2,5)
  plot(ACC,acc(2,:),'--','color',col(1,:),'linewidth',2.5)
  hold on
  plot(ACC,acc(1,:),'color',col(3,:),'linewidth',2.5)
  hold on
  plot(ACC,acc(3,:),':','color',col(4,:),'linewidth',2.5)
  hold on
  plot(ACC,acc(7,:),'color',col(5,:),'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 2 0 8])
set(gca,'XTick',[0 0.5 1 1.5 2])
set(gca,'yTick',[0 2 4 6 8])
grid on
xlabel('acc. (g)')
title('PFA1')

subplot(3,2,6)
 plot(ACC,acc(5,:),'--','color',col(1,:),'linewidth',2.5)
  hold on
  plot(ACC,acc(4,:),'color',col(3,:),'linewidth',2.5)
  hold on
  plot(ACC,acc(6,:),':','color',col(4,:),'linewidth',2.5)
  hold on
  plot(ACC,acc(8,:),'color',col(5,:),'linewidth',2.5)
  set(gca, 'FontName', 'Times')
set(gca, 'FontSize', 15)
axis([0 2 0 8])
set(gca,'XTick',[0 0.5 1 1.5 2])
set(gca,'yTick',[0 2 4 6 8])
grid on
xlabel('acc. (g)')
title('PFA4')
legend('PGA','SaT1','SaT2','SaT3')