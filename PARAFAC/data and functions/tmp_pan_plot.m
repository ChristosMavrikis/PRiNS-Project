%PARAFAC - first use process_pan then this
X3=squeeze(outp(:,:,3));
X8=squeeze(outp(:,:,8));

%[sc3, ld3,ss3]=pca(mncn(X3));
%[sc8, ld8,ss8]=pca(mncn(X8));

%biplot([sc3(5:27,1) sc3(5:27,2)],'Scores',[ld3(:,1) ld3(:,2)],'VarLabels',var_labels)
%biplot([sc8(5:27,1) sc8(5:27,2)],'Scores',[ld8(:,1) ld8(:,2)],'VarLabels',var_labels)
%Factorsraw = parafac(outp,2);

%figure;hold on

%subplot(1,3,1);
%scbr=Factorsraw{1};
%plot(scbr(:,1),scbr(:,2),'ro')

%subplot(1,3,2);
%ldattr=Factorsraw{2};
%plot(ldattr(:,1),ldattr(:,2),'bo')

%subplot(1,3,3);
%scass=Factorsraw{3};
%plot(scass(:,1),scass(:,2),'go')


%X3wn=nprocess(outp,[1 0 0],[0 0 0]);
%[ssX,Corco,It] = pftest(3,X3wn,4);
%Factorsraw = parafac(X3wn,2);
%figure;hold on

%subplot(1,3,1);
%scbr=Factorsraw{1};
%plot(scbr(:,1),scbr(:,2),'ro')
%text(scbr(:,1),scbr(:,2))
%variable_names (:,1) = [];
%variable_names (:,2) = [];
%--------- works ------------
%attr = variable_names;
%Xms2=[mncn(X3); mncn(X8)];

%[scms2, ldms2,ssms2]=pca(Xms2);

%figure;hold on
%plot(scms2(1:30,1),scms2(1:30,2),'ro')
%plot(scms2(31:60,1),scms2(31:60,2),'gs','markerfacecolor','g')
%text(scms2(1:30,1),scms2(1:30,2),num2str(outp(1:30,1,3)))
%text(scms2(31:60,1),scms2(31:60,2),num2str(outp(1:30,1,8)))
%axis equal

nl=10;
%ldms2n=nl*ldms2;
%quiver(zeros(27,1),zeros(27,1),ldms2n(1:27,1),ldms2n(1:27,2),'color','b')
%text(ldms2n(1:27,1),ldms2n(1:27,2),attr)
%title('multiset, centreren per assessor')
%--------- works

X3wn=nprocess(outp,[1 0 1],[0 0 0]);
%[ssX,Corco,It] = pftest(3,X3wn,4);
Factorsraw = parafac(X3wn,2);
figure;hold on
scbr=Factorsraw{1};
ldattr=nl*Factorsraw{2};
scass=nl*Factorsraw{3};

plot(scbr(:,1),scbr(:,2),'ro')
quiver(zeros(23,1),zeros(23,1),ldattr(5:27,1),ldattr(5:27,2),'bo-')
plot(scass(:,1),scass(:,2),'g<','markerfacecolor','g')

legend({'Product Number';'Attributes';'Panelists'})


text(scbr(:,1),scbr(:,2),num2str(outp(1:30,1,1)))
text(ldattr(5:27,1),ldattr(5:27,2),var_labels)
text(scass(:,1),scass(:,2),num2str([1:8]'))

title('Triplot, Product number, attributes and panelist')
