function [W,pval] =  ttest(X,Y)
%=========================================================================
% Feature ranking according to the Ttest statistic
%=========================================================================   
% A=ttest(X,Y) returns a Ttest object initialized with hyperparameters H. 
%
% The train method ranks features with the Ttest statistic.
% The top ranking features are selected and the new data matrix returned.
% 
% Arguments:
%
%   X           -- Input matrix (m x N): m:samples, N: dimension
%   Y           -- Output vector (m x 1): m: saples
%
%  If several thresholds are provided, all criteria are satisfied,
%  e.g. feature_number <= f_max and W > theta. The maximum number of
%  features is never exceeded.
%
%  Model
%
%  W             -- Ranking criterion abs(T statistic), the larger, the better.  
%                     These values are unsorted. All the values are
%                     returned, not just the ones matching the threshold
%                     criteria.
%  pval          -- Pvalues . These values are unsorted.
%

classes = unique(Y);
if length(classes)>2
    error('Only accept 2 classes')
end

Posidx=find(Y==classes(1));
Negidx=find(Y==classes(2));
Mu1=mean(X(Posidx,:));
Mu2=mean(X(Negidx,:));
n1=length(Posidx);
n2=length(Negidx);

Var1=var(X(Posidx,:));
Var2=var(X(Negidx,:));

% Two different ways of computing the standard error
% depending on whether the classes have the same variance or not
Dfree=n1+n2-2;
Spooled=sqrt(((n1-1)*Var1 + (n2-1)*Var2) / Dfree);
Stderr=Spooled * sqrt(1/n1 + 1/n2);

% The t statistic
Stderr(find(Stderr==0))=1;
Tstat=abs(Mu1-Mu2)./Stderr; % We take the absolute value because it does no matter which site is largest

W=Tstat;
[~,alg.fidx]=sort(-Tstat);
pval = 2*tcdf(-Tstat,Dfree); % One tailed test; The maximum pvalue will be 0.5


end