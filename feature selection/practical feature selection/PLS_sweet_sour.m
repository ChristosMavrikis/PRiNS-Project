
%load data
X = A1(:,6:26);
Y = A1(:,27:28); 

Xvars = variable_names(7:27);
Yvars = ["M-Sour";"M-Sweet"];


[S,RMSECV,RMSEP,model]= plsdaforward(X,Y,10);

[minRMSECV, ind_minRMSECV] = min(RMSECV);
disp('------------')
disp('Cream Cheese dataset')
disp('PLS2-R sequential forward search')
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
title('Cream Cheese dataset - PLS-R sequential forward search')
xticks([1:length(Xvars)]);  xticklabels(Xvars(S))
xtickangle(45)

model.validation = plsvalidate(model);


figure;
for i=1:length(Yvars)
    subplot(1, length(Yvars), i)
    hold on;
    plot(model.calibration.Yr(:, i), model.calibration.Yp(:, i), 'b.', 'MarkerSize', 15);
    plot(xlim, ylim, 'k:');
    text(min(xlim), max(ylim), {[' R^2 = ' num2str(round(corr(Y(:, i), model.calibration.Yp(:, i)).^2, 3))]; [' RMSE = ' num2str(round(sqrt(mean((Y(:, i) - model.calibration.Yp(:, i)).^2)), 3))]}, 'FontWeight', 'bold', 'Color', 'b', 'VerticalAlignment', 'top')
    xlabel('Reference');
    ylabel('Prediction');
    title([Yvars{i} ' (PLS2, Calibration)']);
end


figure;
for i=1:length(Yvars)
    subplot(1, length(Yvars), i)
    hold on;
    plot(Y(:, i), model.validation.Yppred(:, i), 'r.', 'MarkerSize', 15);
    plot(xlim, ylim, 'k:');
    text(min(xlim), max(ylim), {[' R^2 = ' num2str(round(corr(Y(:, i), model.validation.Yppred(:, i)).^2, 3))]; [' RMSE = ' num2str(round(sqrt(mean((Y(:, i) - model.validation.Yppred(:, i)).^2)), 3))]}, 'FontWeight', 'bold', 'Color', 'r', 'VerticalAlignment', 'top')
    xlabel('Reference');
    ylabel('Prediction');
    title([Yvars{i} ' (PLS2, Validated)']);
end
