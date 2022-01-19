function vipScore =  vipscore(pls_model)
%=========================================================================
% Compute the VIP score from the plsmodel
%=========================================================================   
% Compute the VIP score of the plsmodel trained with function plstrain
%
% The method returns the VIP scores.
% 
% Arguments:
%
%   plsmodel    -- pls model trained with function plstrain
%
% Returns
%  
%   vipscore    -- VIP score (matrix N x C); N:number of variables
%                  and C number of classes
%

% Compute vip score
for i=1:size(pls_model.Y,2)
    W0 = pls_model.W ./ sqrt(sum(pls_model.W.^2,1));
    B = pls_model.Q;
    T = pls_model.T;
    SS = B(i,:)'.^2.*diag(T'*T);
    tmp = repmat(SS',size(pls_model.X,2),1);
    vipScore(:,i) = sqrt(size(pls_model.X,2)*sum(W0.^2.*tmp,2)/sum(SS));
end
end