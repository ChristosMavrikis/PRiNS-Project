function [scores,loads,ssq,res,q,tsq,tsqs] = pca(data,plots,scl,lv)
%PCA Principal components analysis
%  PCA uses svd to perform pca on a data matrix. It is
%  assumed that samples are rows and variables are columns. 
%  The input is the data matrix (data). Outputs are the scores
%  (scores), loadings (loads), variance info (ssq), residuals
%  (res), Q limit (reslm), T^2 limit (tsqlm), and T^2's (tsq).
%
%  Optional inputs are (plots) plots = 0 suppresses all plots,
%  plots =  1 [default] produces plots with no confidence limits,
%  plots =  2 produces plots with limits,
%  plots = -1 plots the eigenvalues only (without limits),
%  a vector (scl) for plotting scores against, (if scl = 0 sample 
%  numbers will be used), and a scalar (lv) which specifies the
%  number of principal components to use in the model and
%  which suppresses the prompt for number of PCs.
%
%I/O: [scores,loads,ssq,res,reslm,tsqlm,tsq] = pca(data,plots,scl,lvs);
%
%  Note: with plots = 0 and lv specified, this routine requires
%  no interactive user input. If you would like to scale the data
%  before processing use the functions AUTO or SCALE.
%
%See also: EVOLVFA, EWFA, BIGPCA, PCAGUI, PCAPRO, PLTLOADS, PLTSCRS,
%          SCRPLTR, SIMCA, RESMTX, TSQMTX, PCADEMO

%Copyright Eigenvector Research, Inc. 1991-2000
%BMW: 11/93,12/94,2/95,5/95,8/97,12/99
%NBG: 2/96,3/96,10/96,11/98,3/99
%NBG: 10/27/00 NaN/inf check

if nargin < 2
  plots = 1;
end
if (plots > 2 | plots < -1)
  error('Plot option must be -1, 0, 1, or 2')
end
mlimt = 501;
[m,n] = size(data);
if ~all(all(isfinite(data)))    %nbg 10/27/00, start
  error('Input (data) can not contain ''NaN'' or ''inf''')
end                              %nbg 10/27/00, end
if nargin == 4  
  if (lv > min([m n])) %& (plots ~= 0)
    disp('Resetting lv to be equal to min([m n])')
    lv = min([m n]);
  end
end
if n < m
  cov = (data'*data)/(m-1);
  [u,s,v] = svd(cov);
  temp2 = (1:n)';
  escl = 1:n;
else
  cov = (data*data')/(m-1);
  [u,s,v] = svd(cov);
  v = data'*v;
  for i = 1:m
    v(:,i) = v(:,i)/norm(v(:,i));
  end
  temp2 = (1:m)';
  escl  = 1:m;
end
mns    = mean(data);
ssqmns = mns*mns';
ssqtot = sum(diag(cov));
if (ssqtot/ssqmns<1e10)&(plots~=0)
  disp('  ')
  disp('Warning: Data does not appear to be mean centered.')
  disp(' Variance captured table should be read as sum of')
  disp(' squares captured.') 
end
temp = diag(s)*100/(sum(diag(s)));
ssq  = [temp2 diag(s) temp cumsum(temp)];

%  This section calculates the standard errors of the
%  eigenvalues and plots them
mescl    = length(escl);
if plots == 2
  eigmax = ssq(:,2)/(1-(1.96*sqrt(2/m)));
  eigmin = ssq(:,2)/(1+(1.96*sqrt(2/m))); clf
  if mescl<21
    plot(escl,ssq(:,2),'-r',escl,eigmax,'--b',escl,eigmin,'--b',...
      escl,ssq(:,2),'og')
	  disp('first')
  else
    plot(escl(1:20),ssq(1:20,2),'-r',escl(1:20),eigmax(1:20),'--b',...
	  escl(1:20),eigmin(1:20),'--b',escl(1:20),ssq(1:20,2),'og')
	  disp('second')
  end
  title('Eigenvalue vs. PC Number with 95% Confidence Limits')
  xlabel('PC Number'), ylabel('Eigenvalue')
elseif (plots == 1 | plots == -1)
  clf
  if mescl<21
    plot(escl,ssq(:,2),'-r',escl,ssq(:,2),'og')
  else
    plot(escl(1:20),ssq(1:20,2),'-r',escl(1:20),ssq(1:20,2),'og')
  end
  title('Eigenvalue vs. PC Number')
  xlabel('PC Number'), ylabel('Eigenvalue')
end 
%  Print variance information
if plots~=0|nargin<4
  disp('   ')
  disp('        Percent Variance Captured by PCA Model')
  disp('  ')
  disp('Principal     Eigenvalue     % Variance     % Variance')
  disp('Component         of          Captured       Captured')
  disp(' Number         Cov(X)        This  PC        Total')
  disp('---------     ----------     ----------     ----------')
  format = '   %3.0f         %3.2e        %6.2f         %6.2f';
  mprint = min([20 n m]);
  for i = 1:mprint
    tab = sprintf(format,ssq(i,:)); disp(tab)
  end
  if min([n m]) > mprint
    if (nargin < 4 & plots ~= -1)
      sf = input('Print remaining variance information? (yes = 1)');
      if (~isempty(sf))&(sf == 1)
        for i = mprint+1:min([n m])
          tab = sprintf(format,ssq(i,:)); disp(tab)
        end
      end
    end
  end
end
if nargin < 4
  flag = 0;
  while flag == 0;
    lv = input('How many principal components do you want to keep?  ');
    if lv > min([n m])
      disp('Number of PCs must be <= min of number of samples and variables')
    elseif lv < 1
      disp('Number of PCs must be > 0')
    elseif isempty(lv)
      disp('Number of PCs must be > 0')
    else
      flag = 1;
    end
  end
elseif plots~=0
  sf = sprintf('Now calculating statistics based on %g PC model',lv);
  disp(sf)
end
%  Form the PCA Model Based on the Number of PCs Chosen
loads  = v(:,1:lv);
scores = data*loads;
%  Calculate the standard error on the PC loadings if needed
if plots == 2
  loaderr = zeros(n,lv);
  if n > m, nn = m; else, nn = n; end
  for i = 1:lv
    dif = (ssq(:,2)-ones(nn,1)*ssq(i,2)).^2;
    dif = sort(dif);
    sig = sum((ones(nn-1,1)*ssq(i,2))./dif(2:nn,1));
    loaderr(:,i) = ((ssq(i,2)/m)*loads(:,i).^2)*sig;
  end
  loadmax = loads+loaderr;
  loadmin = loads-loaderr;
end
%  Calculate the residuals matrix and Q values
if lv < min([m n])
  resmat = (data - scores*loads')';
  res    = (sum(resmat.^2))';
else 
  res = zeros(m,1);
  if plots~=0
    disp('Residuals not calculated when number of PCs')
    disp(' Equals the number of samples or variables')
    disp(' (res and q = 0)')
  end
end
