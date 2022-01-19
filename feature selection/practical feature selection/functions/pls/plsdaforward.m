function [S, RMSECV, RMSEP,model] = plsdaforward(X, YC, K)
% Runs the sequential forward search with the PLS-DA
% 
% X: input matrix
% YC: dummy input variable
% K: maximum number of features to select
%
% Outputs:
% Top ranked features 
% RMSECV: cross validation rmse
% RMSEP: predicted rmse

A = 1:size(X,2);
S = [];

for j=1:K
    clear tmp_RMSECV tmp_RMSEP
    for i=1:length(A)
        %model = plstrain(X(:,[S A(i)]), YC, 1, 1, -13, 'optscheme', 'vblinds', 'k', '5');
        model = plstrain(X(:,[S A(i)]), YC, 2, 1, -10, 'optscheme', 'vblinds', 'optk', 10, 'type', 'regression');
        model.validation = plsvalidate(model);
        tmp_RMSECV(i) = model.calibration.RMSE;
        tmp_RMSEP(i) = model.validation.RMSE;
    end
    [min_rmse, ind_rmse] = min(tmp_RMSECV);
    S = [S A(ind_rmse)];
    A(ind_rmse)=[];
    RMSECV(j) = min_rmse;
    RMSEP(j) = tmp_RMSEP(ind_rmse);
end
end