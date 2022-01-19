function [Pr] = pca_cv(X,pc,scl,Lo,Lv, hussle)
% [Pr] = pca_cv(X,pc,scl,Lo,Lv, hussle)
%
%  Purpose:
%  PCA Cross validation by leaving out columns and rows
% (Eastment and Krzanowski Technometrics 1982 & Louwerse,
%  Kiers and Smilde Technometrics note 1998).
%
%  Parameters:
%  X        -   Data matrix (n_objecten x n_variables).
%  pc       -   (optional) Maximum number of PC to be used.
%               Default: min( (number of vars - Lv),(number of objects - Lo))
%  scl      -   (optional) Scaling used for data matrix X.
%               0: no scaling
%               1: column centering
%               2: column auto scaling
%               Default: 0
%  Lo       -   (optional) Number of objects that will be left out.
%               Default: 1
%  Lv       -   (optional) Number of variables that will be left out.
%               Default: 1
%  hussle   -   (optional) Randomize variables before cross validation when Lv > 1.
%               1: randomize
%               0: do not randomize
%               Default: 0
%
%  Returned parameters:
%  Pr       -   PRESS as a function of pc number.
%         
%

%   August 1997, University of Amsterdam
%
%   Ad Louwerse
%   Laboratory for Analytical Chemistry
%   Process Analysis & Chemometrics
%   Nieuwe Achtergracht 166
%   1018 WV Amsterdam
%
%   Phone: +31 20 525 6547
%   Email: louwerse@anal.chem.uva.nl

% -------------------------- Argument checking ---------------------
if ( nargin < 1 )
    error('Call is: [Pr] = pca_cv(X,pc,scl,Lo,Lv, hussle)');
end

[No,Nv] = size(X);

if (nargin < 6)
    hussle =0;
end
if (hussle ~= 1) & (hussle ~= 0)
    error('hussle parameter should be 1 or 0')
end
if ( nargin < 4 )
    Lo = 1;
end
if ( nargin < 5 )
    Lv = 1;
end
if ( isempty(Lo) )
    Lo=1;
end
if ( isempty(Lv) )
    Lv=1;
end
if ( nargin < 3 )
    scl = 0;
end
if ( nargin < 2 )
    pc = min([No-Lo,Nv-Lv]);
end
if ( isempty(pc) )
    pc = min([No-Lo,Nv-Lv]);
end
if ( pc > (Nv - Lv) )
    error('Number of PCs must be <= (total - left out) number of variables')
end
if ( pc > (No - Lo) )
    error('Number of PCs must be <= (total - left out) number of objects')
end
if round(Nv/Lv) ~= (Nv/Lv)
    error('the number of combinations of left out variable is not an integer')
end
if round(No/Lo) ~= (No/Lo)
    error('the number of combinations of left out objects is not an integer')
end

% -------------------------- End argument checking ---------------------

% If Nv < No use the kernel PCA method with svd otherwise not
% ChemIntLabSys, 36 (1997) 165-172.
%
kern_svd = (Nv - Lv) > (No - Lo);

% Randomize variables
if (hussle == 1) & (Lv > 1)
    mix = randperm(Nv);
    X = X(:,mix);
end

%if scale option is used the data have to be scaled en rescaled
if scl == 1
    X = mncn(X);
end
if scl ==2
    X = auto(X);
end

% singular value correction factors for left out objects and left out variables
Co = sqrt(No/(No - Lo));
Cv = sqrt(Nv/(Nv - Lv));

% number of object slices and variable slices
Nlo = No/Lo;
Nlv = Nv/Lv;

% Calculate the SVD of X to correct sign later on
if kern_svd
    [U0,D0]=eig(X*X');
    [D0,index]=sort(-diag(D0));
    S=diag(sqrt(-D0(1:pc)));
    Uc=U0(:,index(1:pc));
    Vc=X'*Uc/S;
else
    [V0,D0]=eig(X'*X);
    [D0,index]=sort(-diag(D0));
    S=diag(sqrt(-D0(1:pc)));
    Vc=V0(:,index(1:pc));
    Uc=X*Vc/S;
end

% calculate svd models for left out objects
% The square root of S is calculated as two S's are combined later on
disp(['calculating models for left out objects'])
for i = 1:Nlo

    Xl = [X(1:(i-1)*Lo,:);X(i*Lo+1:No,:)];

    if scl == 1
        Xl = mncn(Xl);
    end
    if scl == 2
        Xl = auto(Xl);
    end

    if kern_svd
        [U0,D0]=eig(Xl*Xl');
        [D0,index]=sort(-diag(D0));
        S=diag(sqrt(-D0(1:pc)));
        U=U0(:,index(1:pc));
        V=Xl'*U/S;
    else
        [V0,D0]=eig(Xl'*Xl);
        [D0,index]=sort(-diag(D0));
        S=diag(sqrt(-D0(1:pc)));
        V=V0(:,index(1:pc));
    end
    V_sign = diag(sign(diag(V'*Vc)));
    eval(['V',int2str(i),' = V*V_sign;'])
    eval(['So',int2str(i),' = sqrt(S.*Co);'])
end

% calculate svd models for left out variables
% The square root of S is calculated as two S's are combined later on
disp(['calculating models for left out variables'])
for j = 1:Nlv

    Xl = [X(:,1:(j-1)*Lv) X(:,j*Lv+1:Nv)];

    if kern_svd
        [U0,D0]=eig(Xl*Xl');
        [D0,index]=sort(-diag(D0));
        S=diag(sqrt(-D0(1:pc)));
        U=U0(:,index(1:pc));
    else
        [V0,D0]=eig(Xl'*Xl);
        [D0,index]=sort(-diag(D0));
        S=diag(sqrt(-D0(1:pc)));
        V=V0(:,index(1:pc));
        U=Xl*V/S;
    end
    U_sign = diag(sign(diag(U'*Uc)));
    eval(['U',int2str(j),' = U*U_sign;'])
    eval(['Sv',int2str(j),' = sqrt(S.*Cv);'])
end

disp(['calculating PRESS'])
Pr=zeros(1,pc);
for i = 1:Nlo
    Is = (i-1)*Lo+1;
    Ie   = i*Lo;
    for j = 1:Nlv
        Js = (j-1)*Lv+1;
        Je   = j*Lv;
        for k = 1:pc
%           calulate: Xe = U[j](:,1:k) * Sv[j](1:k,1:k) * So[i](1:k,1:k) * V[i](:,1:k)';
            eval(['Xe = U',int2str(j),'(Is:Ie,1:k) * Sv',...
                 int2str(j),'(1:k,1:k) * So',int2str(i),'(1:k,1:k) * V',...
                 int2str(i),'(Js:Je,1:k)'';'])
            Xdiff = Xe - X(Is:Ie,Js:Je);
            Pr(k) = Pr(k) + sum(sum(Xdiff.^2));
        end
    end
end
