%Factorsraw = parafac(outp(:,6:27,:),2);
tmp_names  = variable_names;

tmp_names(:,1) = [];
tmp_names(:,2) = [];

names = ["16%" ;"16%" ;"16%";"24%";"24%";"24%";"33%";"33%";"33%";"A-Prot";"A-Prot";"A-Prot";"A-Prot";"A-Prot";"A-Prot";
   "B-Prot";"B-Prot";"B-Prot";"C-Cho";"C-Cho";"C-Cho";"D-Cho";"D-Cho";"D-Cho";"P";"P";"P";"P+Aroma";"P+Aroma";"P+Aroma"];
%figure;hold on

%subplot(1,3,1);
%scbr=Factorsraw{1};
%plot(scbr(:,1),scbr(:,2),'ro')
%text(scbr(:,1),scbr(:,2),names)
%title('Product Number')

%subplot(1,3,2);
%ldattr=Factorsraw{2};
%plot(ldattr(:,1),ldattr(:,2),'bo')
%text(ldattr(:,1),ldattr(:,2),tmp_names(6:27))
%title('attributes')

%subplot(1,3,3);
%scass=Factorsraw{3};
%plot(scass(:,1),scass(:,2),'go')
%text(scass(:,1),scass(:,2),num2str([1:8]'))
%title('Panelists')
%Suptitle('parafac raw data')

%-----------------------
X3wn=nprocess(outp(:,5:27,:),[1 0 0],[0 0 0]);
[ssX,Corco,It] = pftest(3,X3wn,4);
Factorsraw = parafac(X3wn,2);
figure;hold on

subplot(1,3,1);
scbr=Factorsraw{1};
plot(scbr(:,1),scbr(:,2),'ro')
text(scbr(:,1),scbr(:,2),names)
title('Product Name')

subplot(1,3,2);
ldattr=Factorsraw{2};
plot(ldattr(:,1),ldattr(:,2),'bo')
text(ldattr(:,1),ldattr(:,2),tmp_names(5:27))
title('Sensory Variables')

subplot(1,3,3);
scass=Factorsraw{3};
plot(scass(:,1),scass(:,2),'go')
text(scass(:,1),scass(:,2),num2str([1:8]'))
title('Panelists')
Suptitle('parafac data centered across product number')