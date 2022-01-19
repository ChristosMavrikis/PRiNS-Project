%process that data so it can be used by PARAFAC functions from tutorial.
%Run this before running tmp_pan_plot - which plots a triplot from tutorial
 pan1 = find(1==A1(1:240,2));
 pan2 = find(2==A1(1:240,2));
 pan3 = find(3==A1(1:240,2));
 pan4 = find(4==A1(1:240,2));
 pan5 = find(5==A1(1:240,2));
 pan6 = find(6==A1(1:240,2));
 pan7 = find(7==A1(1:240,2));
 pan8 = find(8==A1(1:240,2));

 a_pan1 = A1(pan1,:);
 a_pan2 = A1(pan2,:);
 a_pan3 = A1(pan3,:);
 a_pan4 = A1(pan4,:);
 a_pan5 = A1(pan5,:);
 a_pan6 = A1(pan6,:);
 a_pan7 = A1(pan7,:);
 a_pan8 = A1(pan8,:);

 a_pan1(:,2) = [];
 a_pan2(:,2) = [];
 a_pan3(:,2) = [];
 a_pan4(:,2) = [];
 a_pan5(:,2) = [];
 a_pan6(:,2) = [];
 a_pan7(:,2) = [];
 a_pan8(:,2) = [];
 outp= zeros(30, 27, 8);
 outp(:,:,1) = a_pan1;
 outp(:,:,2) = a_pan2;
 outp(:,:,3) = a_pan3;
 outp(:,:,4) = a_pan4;
 outp(:,:,5) = a_pan5;
 outp(:,:,6) = a_pan6;
 outp(:,:,7) = a_pan7;
 outp(:,:,8) = a_pan8;