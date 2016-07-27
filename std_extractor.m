function [std_vec,number] = std_extractor(data)
dummy = size(data);
rows = dummy(1);
cols = dummy(2);
for ii = 1:cols
K = data(:,ii);
K = K(K~=0);
K = K(K<=0.08);
number(ii) = max(size(K));
std_vec(ii) = (log(std((K))));
end
end