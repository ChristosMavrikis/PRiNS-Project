X3=squeeze(outp(:,16:27,4));
X8=squeeze(outp(:,16:27,8));
tmp_names  = variable_names;

tmp_names(:,1) = [];
tmp_names(:,2) = [];

Xmb=[X3, X8];
Xms=[X3; X8];

[scmb, ldmb,ssmb]=pca(mncn(Xmb),0,[],2);
[scms, ldms,ssms]=pca(mncn(Xms),0,[],2);
names = ["16%" ;"16%" ;"16%";"24%";"24%";"24%";"33%";"33%";"33%";"A-Prot";"A-Prot";"A-Prot";"A-Prot";"A-Prot";"A-Prot";
   "B-Prot";"B-Prot";"B-Prot";"C-Cho";"C-Cho";"C-Cho";"D-Cho";"D-Cho";"D-Cho";"P";"P";"P";"P+Aroma";"P+Aroma";"P+Aroma"];

figure;hold on
plot(scmb(:,1),scmb(:,2),'ro')
text(scmb(:,1),scmb(:,2),names)
axis equal

nl=10;
ldmbn=nl*ldmb;

quiver(zeros(12,1),zeros(12,1),ldmbn(1:12,1),ldmbn(1:12,2),'color','r')
quiver(zeros(12,1),zeros(12,1),ldmbn(13:24,1),ldmbn(13:24,2),'color','g')
text(ldmbn(1:12,1),ldmbn(1:12,2),tmp_names(:,16:27))
text(ldmbn(13:24,1),ldmbn(13:24,2),tmp_names(:,16:27))
title('Multiblock');
legend({'Product number';'Panelist 4';'Panelist 8'});