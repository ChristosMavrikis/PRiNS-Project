function [S,tElapsed]=mRMRGaussian(X,y,K)
%=========================================================================
% Select top K features from the mRMR criteria using Gaussian distribution
%=========================================================================   
%
% The method returns the VIP scores.
% 
% Arguments:
%
%   X    -- input data (m x N)
%   y    -- output data (m x 1)
%   K    -- number of features to select
%
% Returns
%  
%   S	 -- Selected features
%

tStart=tic;

A=1:size(X,2); % Index input features
S=[]; % Index selected features

for i=1:size(X,2)
    rel(i)=miG(X(:,i),y);
end
[~,indmax] = max(rel);
S(1)=A(indmax);
A(indmax) = [] ;


for h=1:K
    diff=zeros(1,length(A));
    for i=1:length(A)
        red=0;
        for j=1:length(S)
            red(j)=miG(X(:,A(i)),X(:,S(j)));
        end
        diff(i) = rel(A(i))-mean(red);
    end

    [~,indmax] = max(diff);
    S(h)=A(indmax);
    A(indmax) = [] ;
end

tElapsed=toc(tStart);

end

function mi = miG(X,Y)
    mi = -1/2*log(1-corr(X,Y)^2);
end


