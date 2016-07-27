function [IM,EDP] = data_selection(data_array,number)
Sa = data_array(:,1);
Dr = data_array(:,2);
bins = 0.05:0.05:3;
IM = zeros(1,1);
EDP = zeros(1,1);
for ii = 1:max(size(bins))-1
lower = bins(ii);
upper = bins(ii+1);
req = find(Sa>lower & Sa<upper);
if max(size(req))>=number
    random_index = randsample(max(size(req)),number);
    K = req(random_index);
    K1 = Sa(K);
    K2 = Dr(K);
else
    K1 = Sa(req);
    K2 = Dr(req);
end
IM = vertcat(IM,K1);
EDP = vertcat(EDP,K2);
end
IM = IM(2:max(size(IM)));
EDP = EDP(2:max(size(EDP)));
end