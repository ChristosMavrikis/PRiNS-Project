function [Yp, RMSE] = plsforward(X, Y, K)
% Runs the sequential forward search 

A = 1:size(X,2);
S = [];
K = 13;
for j=1:K
    clear tmp_RMSECV tmp_RMSEP
    for i=1:length(A)
        model = plstrain(X(:,[S A(i)]), YC, 1, 1, -10, 'optscheme', 'vblinds', 'k', '5');
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