function biplotx(scores, loadings, ssq, Data)

figure;hold on
plot(scores(Data.Treatments==1,1),scores(Data.Treatments==1,2),'ro')
plot(scores(Data.Treatments==2,1),scores(Data.Treatments==2,2),'bs')
plot(scores(Data.Treatments==3,1),scores(Data.Treatments==3,2),'g<')

scoremax=max(sqrt(sum(scores.^2,2)));
loadingmax=max(sqrt(sum(loadings.^2,2)));
loadingscale=loadings./loadingmax.*scoremax;

quiver(zeros(size(loadings,1),1),zeros(size(loadings,1),1),loadingscale(:,1),loadingscale(:,2),'k')
text(loadingscale(:,1),loadingscale(:,2), Data.Variables)
eval(['xlabel(''PC 1 (' num2str(round(ssq(1,3))) '%)'')'])
eval(['ylabel(''PC 2 (' num2str(round(ssq(2,3))) '%)'')'])
legend(Data.Treatment_Label)