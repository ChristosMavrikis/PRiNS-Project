% Makes a 3D plot of the profiles in the data set. By default, column
% profiles (variables) are plotted. To plot the row profiles (observations)
% write: plot_data3D(data,'row'). To define the X axis, write 'X' followed
% by a vector.

function plot_data3D(data,varargin)

figure
[m,n] = size(data);

if any((strcmp(varargin,'X')))
    Xaxis = varargin{find(strcmp(varargin,'X'))+1};
else
    Xaxis = 1:n;
end
hold on

if nargin == 1 || strcmp(varargin,'column')
    for i=1:n
        plot3(1:m,i*ones(1,m),data(:,i))
    end
elseif any(strcmp(varargin,'row'))
    for i=1:m
        plot3(Xaxis,i*ones(1,n),data(i,:))
    end
else
    error('Please specify the plotting direction (row or column)')
end
xlabel('observations')
ylabel('variables')
zlabel('intensity')
view(gca,[31.5 43.2]);