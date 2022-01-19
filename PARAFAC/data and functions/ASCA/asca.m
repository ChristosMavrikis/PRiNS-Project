function [Model]=asca(X,DesignMat2,DesignMat3,FacI,FacII,FacIII)

%   Anova - SCA
%   This is the ASCA algorithm for unbalanced data, for balanced data this
%   algorithm generalizes to the algorithm for balanced data.
%
%   inputs:
%       X           ((I x K) x J)     =   data
%       DesignMat   ((I x K) x I)     =   Matrix containing 0's and 1's defining which samples belong to which individual
%       DesignMat2  ((I x K) x C)     =   Matrix containing 0's and 1's defining which samples belong to which class
%       DesignMat3  ((I x K) x K)     =   Matrix containing 0's and 1's defining which samples belong to which measurement occasion
%       FacI                          =   number of factors for submodel K
%       FacII                         =   number of factors for submodel Kh
%       FacIII                        =   number of factors for submodel Khi
%       Sc                            = ` Optional parameter to determine scaling (1=scaling, 0=no scaling, default is no scaling)
%
%   outputs:
%       Model                         =   Structure Array containing the model
%       Model.dimensions              =   Array containing all dimensions of the data and the model
%
%       Model.offset        (1 x J)             = overall offset of the data
%
%       Model.total.expvartot                   =   Total percentage of explained variation 
%
%       Model.I.data                      =   Data used for submodel K
%       Model.I.scores                    =   submodel K scores
%       Model.I.loadings                  =   submodel K loadings
%       Model.I.percexp.total             =   Total percentage of explained variation for submodel K
%       Model.I.percexp.perpc             =   Percentage of explained variance for each PC
% 
%       Model.II.data                      =   Data used for the submodel Kh
%       Model.II.scores                    =   submodel Kh scores
%       Model.II.{class}                   =   Cell Array containing information per treatment group
%       Model.II.loadings                  =   submodel Kh loadings
%       Model.II.percexp.total             =   Total percentage of explained variance
%       Model.II.percexp.perpc             =   Percentage of explained variance for each PC
%
%       Model.III.data                       =   Data used for submodel Khi
%       Model.III.scores                     =   submodel Khi scores
%       Model.III.{ind}                      =   Cell Array containing information per individual
%       Model.III.loadings                   =   Within-individual loadings
%       Model.III.percexp.total              =   Total amount of explained variance 
%       Model.III.percexp.perpc              =   Percentage of explained variance of each PC
%
%       I/O: [Model]=asca(X,DesignMat,DesignMat2,DesignMat3,FacI,FacII,FacIII)


%I   =   sum(DesignMat1);
J   =   size(X,2);
C   =   size(DesignMat2,2);
K   =   size(DesignMat3,2);
Kall=   size(X,1);

%   Part I: Calculate Overall Mean
Model.mean  =   mean(X);

%   Subtract offset from Data
Xmean    =   X-ones(Kall,1)*Model.mean;
Model.X=Xmean;
%   Submodel I
Xclass=[];
for k=1:K
    [smpsK,dummy]=find(DesignMat3(:,k));
    XK(k,:)=mean(Xmean(smpsK,:));
end
Model.I.X=DesignMat3*XK;

[u s v]=svd(Model.I.X);
Model.I.u=u(:,1:FacI);
Model.I.s=s(1:FacI,1:FacI);
Model.I.v=v(:,1:FacI);
Model.I.Scores=Model.I.u*Model.I.s;
Model.I.Loadings=Model.I.v;
Model.I.SSQ=(diag(s).^2)/(sum(diag(s).^2))*100;
clear u s v XK NK smpsK dummy
%   Deflate the data for the II Submodel
Xwitclass=Xmean-Model.I.X;
%   Generate new Design matrix combining DesignMat2 and DesignMat3 
icol=0;
for i2=1:size(DesignMat2,2)
    for i3=1:size(DesignMat3,2)
        icol=icol+1;
        ismpnow=intersect(find(DesignMat2(:,i2)),find(DesignMat3(:,i3)));
        DesignMat4(ismpnow,icol)=1;
        XKh(icol,:)=mean(Xwitclass(ismpnow,:));
    end
end
Model.II.X=DesignMat4*XKh;


[u s v]=svd(Model.II.X);
Model.II.u=u(:,1:FacII);
Model.II.s=s(1:FacII,1:FacII);
Model.II.v=v(:,1:FacII);
Model.II.Scores=Model.II.u*Model.II.s;
Model.II.Loadings=Model.II.v;
Model.II.SSQ=(diag(s).^2)/(sum(diag(s).^2))*100;
clear u s v XKh icol

%   Deflate the data for the III Submodel
Model.III.X=Xwitclass-Model.II.X;
[u s v]=svd(Model.III.X);
Model.III.u=u(:,1:FacIII);
Model.III.s=s(1:FacIII,1:FacIII);
Model.III.v=v(:,1:FacIII);
Model.III.Scores=Model.III.u*Model.III.s;
Model.III.Loadings=Model.III.v;
Model.III.SSQ=(diag(s).^2)/(sum(diag(s).^2))*100;

%%  Additional diagnostics
% Total percentage of explained variation
RES =   Xmean-Model.I.Scores*Model.I.Loadings'-Model.II.Scores*Model.II.Loadings'-Model.III.Scores*Model.III.Loadings';
Model.SSQ   =   (1-(ssq(RES)/ssq(Xmean)))*100;

%  Calculate % per class
for c=1:C
    ismp=find(DesignMat2(:,c));
    Model.II.class{c}.X              =   Model.II.X(ismp,:);
    Model.II.class{c}.Scores            =   Model.II.Scores(ismp,:);
    Model.II.class{c}.SSQ.Total     =   (1-(ssq(Model.II.class{c}.X - Model.II.class{c}.Scores*Model.II.Loadings')/ssq(Model.II.class{c}.X)))*100;
    for ifac=1:FacII
        Model.II.class{c}.SSQ.perpc(ifac)     =   (1-(ssq(Model.II.class{c}.X - Model.II.class{c}.Scores(:,ifac)*Model.II.Loadings(:,ifac)')/ssq(Model.II.class{c}.X)))*100;
    end
    
    Model.III.class{c}.X              =   Model.III.X(ismp,:);
    Model.III.class{c}.Scores            =   Model.III.Scores(ismp,:);
    Model.III.class{c}.SSQ.total     =   (1-(ssq(Model.III.class{c}.X - Model.III.class{c}.Scores*Model.III.Loadings')/ssq(Model.III.class{c}.X)))*100;
    for ifac=1:FacIII
        Model.III.class{c}.SSQ.perpc(ifac)     =   (1-(ssq(Model.III.class{c}.X - Model.III.class{c}.Scores(:,ifac)*Model.III.Loadings(:,ifac)')/ssq(Model.III.class{c}.X)))*100;
    end

end

Model.DesignMat=DesignMat4;
Model.DesignMat_time=DesignMat3;
Model.DesignMat_class=DesignMat2;




function t = ssq(a)
%SSQ	SSQ(A) is the sum of squares of the elements of matrix A.
t = sum(sum(a.^2));

function [scores,loads,ssq] = pca(data,lv)
mlimt = 501;
[m,n] = size(data);
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
temp = diag(s)*100/(sum(diag(s)));
ssq  = [temp2 diag(s) temp cumsum(temp)];

%  Form the PCA Model Based on the Number of PCs Chosen
loads  = v(:,1:lv);
scores = data*loads;
