
[mcx,mx] = mncn(A1(:,6:28));
[ax,mx,stdx] = auto(mcx);
[scores,loads,ssq,res] = pca(mcx,1);

plot(loads(:,4),loads(:,3),'ro')
text(loads(:,4),loads(:,3),variable_names(:,7:29))
title('Sensory Variables')
xlabel("Loadings Tree")
ylabel("Loadings Four")

