
X = A1(:,2:28);
[m,N]=size(X);
Xvars = variable_names(:,3:29);
Y = cheese_types();
Yvars = ["16%" "24%" "33%"  "A-Prot" "B-Prot" "C-CHO" "D-CHO" "P" "P+Aroma"];


YCvars = unique(Y, 'stable');
YC = false(size(Y, 1), length(YCvars));
for i=1:length(YCvars)
    YC(:, i) = strcmp(Y, YCvars(i));
end
YC = double(YC);


model_plsda = plstrain(X, YC, 2, 1, -10, 'optscheme', 'vblinds', 'optk', 5, 'type', 'regression');
model_plsda.validation = plsvalidate(model_plsda); 
disp('------------')
disp('Cream Cheese dataset')
disp('PLS-DA - All variables')
disp(['RMSECV: ' num2str(model_plsda.calibration.RMSE)])
disp(['RMSEP: ' num2str(model_plsda.validation.RMSE)])


Yvars = unique(Y);
for i=1:size(YC,2)
    [W,pval] = ttest(X,YC(:,i));
    % plot
    figure(3+i)
    subplot(211)
    bar(W); xticks([1:length(Xvars)]);  xticklabels(Xvars)
    xtickangle(45)
    grid on
    title([Yvars{i} ' vs others'])
    subplot(212)
    bar(pval); xticks([1:length(Xvars)]);  xticklabels(Xvars)
    xtickangle(45)
    grid on
    title('p-value')
end