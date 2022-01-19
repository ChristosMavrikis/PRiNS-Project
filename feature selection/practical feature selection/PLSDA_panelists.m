X = A1(:,6:28);
[m,N]=size(X);
Xvars = variable_names(:,7:29);
Y = A1(:,2);
Yvars = [1 2 3 4 5 6 7 8 ];

YCvars = unique(Y, 'stable');
YC = false(size(Y, 1), length(YCvars));
for i=1:length(YCvars)
    %YC(:, i) = strcmp(Y, YCvars(i));
     YC(:, i) = Y == YCvars(i);
end
YC = double(YC);

% 2.2. Train a PLS-DA model with all variables
model_plsda = plstrain(X, YC, 2, 1, -10, 'optscheme', 'vblinds', 'optk', 5, 'type', 'regression');
model_plsda.validation = plsvalidate(model_plsda); 
disp('------------')
disp('Wine dataset')
disp('PLS-DA - All variables')
disp(['RMSECV: ' num2str(model_plsda.calibration.RMSE)])
disp(['RMSEP: ' num2str(model_plsda.validation.RMSE)])
Yvars = unique(Y);


for i=1:size(YC,2)
    [W,pval] = ttest(X,YC(:,i));
    figure(3+i)
    subplot(211)
    bar(W); xticks([1:length(Xvars)]);  xticklabels(Xvars)
    xtickangle(45)
    grid on
    title([i])
    subplot(212)
    bar(pval); xticks([1:length(Xvars)]);  xticklabels(Xvars)
    xtickangle(45)
    grid on
    title('p-value')
end

figure;
for i=1:length(YCvars)
    subplot(1, length(YCvars), i)
    hold on;
    plot(YC(:, i), model_plsda.validation.Yppred(:, i), 'r.', 'MarkerSize', 15);
    plot(xlim, [0.5 0.5], 'k:')
    xlabel('Reference');
    ylabel('Prediction');
    title(YCvars(i));
end
