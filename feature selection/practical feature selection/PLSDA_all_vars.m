X = A1(:,6:28);
[m,N]=size(X);
Xvars = variable_names(:,7:29);
Y = cheese_types;
Yvars = ["16%" "24%" "33%"  "A-Prot" "B-Prot" "C-CHO" "D-CHO" "P" "P+Aroma"];

YCvars = unique(Y, 'stable');
YC = false(size(Y, 1), length(YCvars));
for i=1:length(YCvars)
    YC(:, i) = strcmp(Y, YCvars(i));
end
YC = double(YC);


model_plsda = plstrain(X, YC, 2, 1, -10, 'optscheme', 'vblinds', 'optk', 5, 'type', 'regression');
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
    [~, j] = max(model_plsda.validation.Yppred(i, :));
    YCpred(i) = YCvars(j);
end

figure;
confusionchart(string(Y), YCpred);
title('Validated PLS-DA results');