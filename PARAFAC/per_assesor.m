X3=squeeze(outp(:,:,3));
X8=squeeze(outp(:,:,8));

[sc3, ld3,ss3]=pca(mncn(X3),0,[],2);
[sc8, ld8,ss8]=pca(mncn(X8),0,[],2);


figure;hold on
plot(sc3(:,1),sc3(:,2),'ro')
%text(sc3(:,1),sc3(:,2),num2str(salt))
axis equal

nl=5;
ld3n=nl*ld3;
quiver(zeros(11,1),zeros(11,1),ld3n(:,1),ld3n(:,2))
%text(ld3n(:,1),ld3n(:,2),attrib)
title('assessor 3')


figure;hold on
plot(sc8(:,1),sc8(:,2),'ro')
text(sc8(:,1),sc8(:,2),num2str(salt))
axis equal

nl=5;
ld8n=nl*ld8;
quiver(zeros(11,1),zeros(11,1),ld8n(:,1),ld8n(:,2))
text(ld8n(:,1),ld8n(:,2),attrib)
title('assessor 8')