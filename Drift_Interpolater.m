function [output] = Drift_Interpolater()
Drift_matrix = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Drift_matrix.txt');
SF = Drift_matrix(1,:);
IM = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Sa.txt');
for ii = 1:max(size(IM))
IM_value = IM(ii);
IM_range = IM_value*SF;
Dr_range = Drift_matrix(ii+1,:);
Dr_range = Dr_range(Dr_range~=0);
K = size(Dr_range);
col = K(2);
IM_range = IM_range(1:col);
output(ii,:) = interp1(IM_range,Dr_range,0.01:0.01:1.5,'linear','extrap');
end
%output(output>cutoff_drift) = NaN;
end