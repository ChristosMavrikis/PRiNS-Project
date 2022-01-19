%DATA CENTERING using 3 different options (mean, median, control)
% 
% [mcx,mx] = mncn(x)
% 
% returns a mean-centered matrix (mcx), obtained by subtracting the column
% means (mx) from the data matrix (x). 
% 
% 
% [mcx,mx] = mncn(x,'median')
% 
% Same as above, but defines mx as column medians instead of means
% 
% 
% [mcx,mx] = mncn(x,'control',T)

% Defines mx as column means of the control group, defined by the vector T.
% The user is asked which value identifies the control group.
%
% Carlo Bertinetto, Radboud University, 19/11/2018

function [mcx,mx] = mncn(x,varargin)

[m,n] = size(x);

if nargin == 1 || any(strcmp(varargin,'mean'))
    mx    = mean(x);
elseif any(strcmp(varargin,'median'))
    mx    = median(x);
elseif any(strcmp(varargin,'control'))
    if nargin < 3
        error('Please include a vector indicating the treatments')
    end
    c = input('Which value corresponds to the control group?');
    T = varargin{find(strcmp(varargin,'control'))+1};
    cntr = T==c;
    mx    = mean(x(cntr,:));
else
    error('Please enter one of the allowed centering methods: mean, median, control')
end
    
mcx   = (x-mx(ones(m,1),:));
