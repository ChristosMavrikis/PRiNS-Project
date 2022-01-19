
X = A1(:,2:28);
[m,N]=size(X);
Xvars = variable_names(:,3:29);
Y = cheese_types(:,1);
Yvars = ["16%" "24%" "33%"  "A-Prot" "B-Prot" "C-CHO" "D-CHO" "P" "P+Aroma"];

YCvars = unique(Y, 'stable');
YC = false(size(Y, 1), length(YCvars));
for i=1:length(YCvars)
    YC(:, i) = strcmp(Y, YCvars(i));
end
YC = double(YC);


[S,RMSECV,RMSEP,model] = plsdaforward(X,YC,18);

[minRMSECV, ind_minRMSECV] = min(RMSECV);
disp('------------')
disp('Cream Cheese dataset')
disp('PLS-DA sequential forward search')
disp(['RMSECV: ' num2str(RMSECV(ind_minRMSECV))])
disp(['RMSEP: ' num2str(RMSEP(ind_minRMSECV))])

figure(7)
plot(RMSECV,'-o','LineWidth',2)
hold on
plot(RMSEP,'-s','LineWidth',2)
grid on
legend('RMSECV','RMSEP')
xlabel('Top ranked features')
ylabel('RMSE')
title('Cream Cheese dataset - PLS-DA sequential forward search')
xticks([1:length(Xvars)]);  xticklabels(Xvars(S))
xtickangle(45)

model.validation = plsvalidate(model);

figure;
for i=1:length(YCvars)
    subplot(1, length(YCvars), i)
    hold on;
    plot(YC(:, i), model.validation.Yppred(:, i), 'r.', 'MarkerSize', 15);
    plot(xlim, [0.5 0.5], 'k:')
    xlabel('Reference');
    ylabel('Prediction');
    title(YCvars(i));
end

predictions = double(model.validation.Yppred>=0.5)


YCpred = strings(size(Y));
for i=1:size(Y, 1)
    [~, j] = max(model.validation.Yppred(i, :));
    YCpred(i) = YCvars(j);
end

figure;
confusionchart(string(Y), YCpred);
title('Validated PLS-DA  Forward results');