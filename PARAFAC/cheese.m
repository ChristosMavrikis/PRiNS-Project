%get the values in the Excel using xlsread.
%Loads data into A1 and excludes  first column and adds 
% variables names into var_labels

close all
clear all
clc
load cheese.mat

%disp(str2double(c(:,1)))

I = load ('cheese.mat');

disp('Data-set has 241 samples & 29 attributes');

variable_names = (I.c(1,:));
tmp_I = I.c;
I.c(1,:) = [];

cheese_types = (I.c(:,1));
I.c(:,1) = [];
for c = 2:28
   fprintf('%s has mean %f and std %f \n',variable_names{c},mean(cellfun(@sum,I.c(1:240,c))),std(cellfun(@sum,I.c(1:240,c))));
end
%cellfun(@plot,I.c(1:240,2));

A1 = I.c(1:240,1:28);
A1 = cell2mat(A1);
disp((A1));
var_labels = variable_names((6:28));
