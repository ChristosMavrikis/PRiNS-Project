X = A1(:,2:28);
[m,N]=size(X);
tmp = variable_names;
tmp (:,1:2) = [];
Xvars = tmp;
Y = cheese_types();
Yvars = ["16%" "24%" "33%"  "A-Prot" "B-Prot" "C-CHO" "D-CHO" "P" "P+Aroma"];
YCvars = unique(Y, 'stable');
YC = false(size(Y, 1), length(YCvars));
for i=1:length(YCvars)
    YC(:, i) = strcmp(Y, YCvars(i));
end
YC = double(YC);


model_plsda = plstrain(X, YC, 2, 1, -13, 'optscheme', 'vblinds', 'optk', 5, 'type', 'regression');
vipScore = vipscore(model_plsda);
indVIP = (vipScore >= 1.15);
indVIP = sum(indVIP,2)>0;

model_plsda_vip = plstrain(X(:,indVIP), YC, 2, 1, -13, 'optscheme', 'vblinds', 'optk', 5, 'type', 'regression');
model_plsda_vip.validation = plsvalidate(model_plsda_vip); 
disp('------------')
disp('Ceram Cheese dataset')
disp('PLS-DA VIP')
disp(['RMSECV: ' num2str(model_plsda_vip.calibration.RMSE)])
disp(['RMSEP: ' num2str(model_plsda_vip.validation.RMSE)])

model_plsda.validation = plsvalidate(model_plsda);
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

YCpred = strings(size(Y));
for i=1:size(Y, 1)
    [~, j] = max(model.validation.Yppred(i, :));
    YCpred(i) = YCvars(j);
end

figure;
confusionchart(string(Y), YCpred);
title('Validated PLS-DA VIP results');