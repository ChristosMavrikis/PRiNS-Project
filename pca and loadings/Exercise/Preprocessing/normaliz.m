function [xn,xsum]=normaliz(x)
% this function normalizes each sample (row) in matrix x to a sum of 1
% matrix xn contains the normalized data
% column vector xsum contains the row sums of the original data x

xsum=sum(abs(x),2);
XSUM=repmat(xsum,1,size(x,2));
xn=x./XSUM;