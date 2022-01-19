%biplot for pca that shows each class
figure()
[mcx,mx] = mncn(A1);
[scores,loads,ssq,res] = pca(mcx(:,2:23),1);
clusters = cheese_types;
% Store handle to biplot
h = biplot([loads(:,3) loads(:,2)],'Scores',[scores(:,3) scores(:,2)],'Varlabels',var_labels(:,2:23));
% Identify each handle
hID = get(h, 'tag'); 
% Isolate handles to scatter points
hPt = h(strcmp(hID,'obsmarker')); 
% Identify cluster groups
grp = findgroups(clusters);    %r2015b or later - leave comment if you need an alternative
grp(isnan(grp)) = max(grp(~isnan(grp)))+1; 
grpID = 1:max(grp); 
class_string = ["16%","24%","33%","A-Prot","B-Prot","C-CHO","D-CHO","P-","P+Aroma"];
% assign colors and legend display name
clrMap = lines(length(unique(grp)));   % using 'lines' colormap
for i = 1:max(grp)
    %set(hPt(grp==i), 'Color', clrMap(i,:),'MarkerSize',14 ,'DisplayName', sprintf('no. Cluster %d Product %s', grpID(i),class_string(i)))
    set(hPt(grp==i), 'Color', clrMap(i,:),'MarkerSize',14 ,'DisplayName', sprintf('Product %s',class_string(i)))
end
% add legend to identify cluster
[~, unqIdx] = unique(grp);
legend(hPt(unqIdx));
title("Cream Cheese Data Mean Centered Biplot")
xlabel("Component 3 10.82% variance")
ylabel("Component 2 13.90% variance")
