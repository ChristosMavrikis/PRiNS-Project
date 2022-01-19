%AUTOSCALING to unit variance and mean (or median or control-mean) zero
% 
% [ax,mx,stdx] = auto(x)
% 
% returns a mean-centered matrix (ax), obtained by subtracting the column
% means (mx) from the data matrix (x), and dividing each column by its
% standard deviation (stdx).
% 
% 
% [ax,mx,stdx] = auto(x,'median')
% 
% Same as above, but defines mx as column medians instead of means.
% 
% 
% [ax,mx,stdx] = auto(x,'control',T)

% Defines mx as column means of the control group, defined by the vector T.
% The user is asked which value identifies the control group.
% 
% 
% [ax,mx,stdx] = auto(x,'omitnan')
% 
% Ignores NaN values during the calculation. Specify the centering method
% in that case.
% 
% 
% Carlo Bertinetto, Radboud University, 19/11/2018


function [ax,mx,stdx] = auto(x,varargin)

[m,n] = size(x);

if any(strcmp(varargin,'omitnan'))
    NaNflag = 'omitnan';
else
    NaNflag = 'includenan';
end

if nargin == 1 || any(strcmp(varargin,'mean'))
    mx    = mean(x,NaNflag);
elseif any(strcmp(varargin,'median'))
    mx    = median(x,NaNflag);
elseif any(strcmp(varargin,'control'))
    if nargin < 3
        error('Please include a vector indicating the treatments')
    end
    c = input('Which value corresponds to the control group?');
    T = varargin{find(strcmp(varargin,'control'))+1};
    cntr = T==c;
    mx    = mean(x(cntr,:),NaNflag);
else    
    error('Please enter one of the allowed centering methods: mean, median, control')
end

stdx  = std(x,NaNflag);
ax    = (x-mx(ones(m,1),:))./stdx(ones(m,1),:);

