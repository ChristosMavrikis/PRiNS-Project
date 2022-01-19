load bread
X3wn=nprocess(X3w,[1 0 1],[0 0 0]);
Factors = parafac(X3wn,2);

figure;hold on
plot(Factors{1}(:,1),Factors{1}(:,2),'ro')
%text(Factors{1}(:,1),Factors{1}(:,2),num2str(salt))

n3=50;
f3=n3*Factors{3};
plot(f3(:,1),f3(:,2),'bx')
text(f3(:,1),f3(:,2),num2str([1:8]'))

n2=50;
f2=n2*Factors{2};
quiver(zeros(11,1),zeros(11,1),f2(:,1),f2(:,2),'g')
text(f2(:,1),f2(:,2),attrib)
