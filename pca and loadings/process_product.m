%process that data so it can be used by PARAFAC functions from tutorial.
%Run this before running tmp_pan_plot - which plots a triplot from tutorial
 prod1 = find(1==A1(1:240,1));
 prod2 = find(2==A1(1:240,1));
 prod3 = find(3==A1(1:240,1));
 prod4 = find(4==A1(1:240,1));
 prod5 = find(5==A1(1:240,1));
 prod6 = find(6==A1(1:240,1));
 prod7 = find(7==A1(1:240,1));
 prod8 = find(8==A1(1:240,1));
 prod9 = find(9==A1(1:240,1));

 a_prod1 = A1(prod1,:);
 a_prod2 = A1(prod2,:);
 a_prod3 = A1(prod3,:);
 a_prod4 = A1(prod4,:);
 a_prod5 = A1(prod5,:);
 a_prod6 = A1(prod6,:);
 a_prod7 = A1(prod7,:);
 a_prod8 = A1(prod8,:);
 a_prod9 = A1(prod9,:);

 a_prod1(:,2) = [];
 a_prod2(:,2) = [];
 a_prod3(:,2) = [];
 a_prod4(:,2) = [];
 a_prod5(:,2) = [];
 a_prod6(:,2) = [];
 a_prod7(:,2) = [];
 a_prod8(:,2) = [];
 a_prod9(:,2) = [];
 
 outp= zeros(24, 27, 8);
 outp(:,:,1) = a_prod1;
 outp(:,:,2) = a_prod2;
 outp(:,:,3) = a_prod3;
 outp(:,:,4) = a_prod4;
 outp(:,:,5) = a_prod5;
 outp(:,:,6) = a_prod6;
 outp(:,:,7) = a_prod7;
 outp(:,:,8) = a_prod8;
 outp(:,:,9) = a_prod9;