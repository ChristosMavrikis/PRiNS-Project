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


[S,RMSECV,RMSEP,model] = plsdaforward(X,YC,10);

[minRMSECV, ind_minRMSECV] = min(RMSECV);
disp('------------')
disp('Wine dataset')
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
title('Wine dataset - PLS-DA sequential forward search')
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