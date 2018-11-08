clear;
close all;

Raw_Data = load('data.parameters.norm.dat');
dvar = load('data.dependent_var.norm.dat'); % Dependent variable
N = size(Raw_Data,1);
d = size(Raw_Data,2);
C = cov(Raw_Data);

% Calculate new base V
[V,D] = eig(C);

eigenvalue = eig(C);

% By default, the "eigenvalue" is in increasing order,
% for visually convience, we want it to be decreasing order.
eigval = zeros(d,1);
for i = 1:d-1
    eigval(i) = eigenvalue(d-i);
end


% Seclect only 3 base vectors
base = V(:,d-2:d);

% Project the data onto the new axis
Proj_Data = transpose(base) * transpose(Raw_Data);
Proj_Data = Proj_Data';

%% Prepare the data
% Randomize data order
msize = length(Proj_Data);
xpermidx = randperm(msize);
data = Proj_Data(xpermidx,:);
dvar = dvar(xpermidx);
x = data(:,1);
y = data(:,2);
z = data(:,3);

%% Initialize the map
% size?
% Through calculation, x[min,max,mean,std] = [-8, 5.2, -1.4, 1.2]
% and y[min,max,mean,std] = [-6.2, 1.8, -3.8, 1.4]
x_lower = -7;
x_uppur = 5;
y_lower = -6;
y_uppur = 2;
step = 0.5;
X0 = x_lower:step:x_uppur;
Y0 = y_lower:step:y_uppur;
Z0 = zeros(length(X0),length(Y0));

% !!!! do not change the boundary of X0 and Y0!!!
% Otherwise the 



%% First plot
% Plot the data
figure;
dp = subplot(3,2,1:4);
for i = 1 : N
    if dvar(i) ==1
        plot3(x(i), y(i), z(i), 'r*');
    elseif dvar(i) ==2
        plot3(x(i), y(i), z(i), 'ro');
    elseif dvar(i) ==3
        plot3(x(i), y(i), z(i), 'c+');
    elseif dvar(i) ==4
        plot3(x(i), y(i), z(i), 'cx');
    elseif dvar(i) ==5
        plot3(x(i), y(i), z(i), 'bs');
    elseif dvar(i) ==6
        plot3(x(i), y(i), z(i), 'kd');
    elseif dvar(i) ==7
        plot3(x(i), y(i), z(i), 'gd');
    end
    hold on;
end

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal;
% Plot map nodes
mesh(X0, Y0, Z0');
row = length(X0);
col = length(Y0);
for i = 1:row
    for j = 1:col
        plot3(X0(i),Y0(j),Z0(i,j),'ok');
    end
end
title('Data Space')
xp = subplot(3,2,5);
title('Learning Parameter')
xlabel ('Iteration Number')

qp = subplot(3,2,6);
title('the value of sigma')
xlabel ('Iteration Number')


%% Start iterations
it = 1;
pt = 0.001; % pause time
%% Initialise learning parameter

Nd = zeros(row*col, 3);
k = 1;
for i = 1:row
    for j = 1:col
        Nd(k, :) = [X0(i) Y0(j) Z0(i, j)]; % the kth row
        k = k +1;
    end
end 

tau = 1;
v = 3;
eta = zeros(6,1);
sigma = zeros(6,1);
eta(1) = exp(-(it/tau));
s0 = 8;
sigma(1) = s0*exp(-(it/v));

cla(xp)
subplot(3,2,5);
bar(eta,'r');
hold on
plot([0 7.5],[0.001 0.001],'k--');
title('Learning Parameter')
xlabel ('Iteration Number')
set(gca,'YScale','log');
set(gca,'YLim',[0.0001 1]); hold off
cla(qp)
subplot(3,2,6);
bar(sigma,'g');
hold on
title('value of sigma')
xlabel ('Iteration Number')
hold off


%% Iterations loop
while eta(it) > 0.001
    %% Find distance of all data points to the nodes
    dists = pdist2(data,Nd);
    %% Initialize data tag
    dataCflag = zeros(msize,1);
    %% loop through data points
    for ii = 1:msize
        dp = subplot(3,2,1:4);
        
        plot3(data(ii,1),data(ii,2),data(ii,3),'dr','markersize',6);
        % Find nearest Node in this point
        [~,p] = min(dists(ii,:));
        plot3(Nd(p,1),Nd(p,2),Nd(p,3),'or','markersize',9);
        %% Move nodes towards this data point
        ve =  data(ii,:)-Nd(p,:);%vector between them
        plot3 ([Nd(p,1),data(ii,1)],[Nd(p,2),data(ii,2)],[Nd(p,3),data(ii,3)],'-g','linewidth',2)
        % plot vector between them
        % Move the nodes
        vect = zeros(length(Nd), 3);
        for jj = 1:length(Nd)
            vect(jj,:) = ve*eta(it)*exp(-((Nd(p,1)-Nd(jj,1))^2 + (Nd(p,2)-Nd(jj,2))^2)/sigma(it)^2); 
            %vector after eta and hfunct aplied
        end 
        quiver3(Nd(p,1),Nd(p,2),Nd(p,3),vect(p,1),vect(p,2),vect(p,3),0,'-r','MaxHeadSize',9,'linewidth',2) ;
        hold off;
        Nd = Nd + vect; %new nodes position
        % update distances for all nodes
        cd = pdist2(data,Nd);
        dists = cd;
        %% Update data cluster tag
        dataCflag(ii) = p; %% p is the index of the node nearst to the data point
        %% Data Iteration Plot update
        % Plot Data
        pause(pt)
        cla(dp)
        dp = subplot(3,2,1:4);
 
        hold on;
        
        
        plot3(data(:,1),data(:,2),data(:,3),'b.','markersize',4); hold on;
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        axis equal;
        % Plot the nodes
        kk = 1;
        X00 = zeros(row, col);
        Y00 = zeros(row, col);
        Z00 = zeros(row, col);
        for i = 1:row
            for j = 1:col
                X00(i,j) = Nd(kk,1);
                Y00(i,j) = Nd(kk,2);
                Z00(i,j) = Nd(kk,3);
                kk = kk +1;
            end
        end
        mesh(X00,Y00,Z00);
        plot3(X00,Y00,Z00,'.r','markersize',9);
        axis equal;
        title('Data Space')
    end
    %% Plot eta and quantasation error
    % Permutate data
    permidx = randperm(msize);
    data = data(permidx,:);
    dvar = dvar(permidx);
    dataCflag = dataCflag(permidx);
    pt = pt*0.1; % increase plot speed for each iteration
    it = it + 1; % increase iteration
     %% Update Learning parameter
    eta(it) = exp(-(it/tau));
    sigma(it) = s0*exp(-it/v);
    cla(qp)
    qp = subplot(3,2,6);
    bar(sigma, 'g');
    title(['value of sigma'])
    xlabel ('Iteration Number')
     pause(pt*20)
    cla(xp)
    subplot(3,2,5);
    bar(eta,'r');hold on
    plot([0 7.5],[0.001 0.001],'k--');
    title('Learning Parameter')
    set(gca,'YScale','log');
    xlabel ('Iteration Number')
    set(gca,'YLim',[0.0001 1]); hold off
end
%% plot the final result
% plot the original data points
figure(2)
%h3 = mesh(X0, Y0, Z0');
hold on;
%set(h3,'facealpha',0.5)
%plot3(data(:,1),data(:,2),data(:,3),'b.','markersize',4);
%hold on;
axis equal;
% plot the final SOM result
%figure(3)
h4 = mesh(X00,Y00,Z00);
hold on;
set(h4,'facealpha',0.5);
%plot3(X00,Y00,Z00,'.r','markersize',9);

for i = 1 : N
    if dvar(i) ==1
        plot3(data(i,1), data(i,2), data(i,3), 'r*');
    elseif dvar(i) ==2
        plot3(data(i,1), data(i,2), data(i,3), 'ro');
    elseif dvar(i) ==3
        plot3(data(i,1), data(i,2), data(i,3), 'c+');
    elseif dvar(i) ==4
        plot3(data(i,1), data(i,2), data(i,3), 'cx');
    elseif dvar(i) ==5
        plot3(data(i,1), data(i,2), data(i,3), 'bs');
    elseif dvar(i) ==6
        plot3(data(i,1), data(i,2), data(i,3), 'kd');
    elseif dvar(i) ==7
        plot3(data(i,1), data(i,2), data(i,3), 'gd');
    end
    hold on;
end
axis equal;

%% Determine the coordinates of the data after projecting on the SOM.
SOM_proj_data = zeros(msize:2);
k = length(Nd);
for ii = 1:msize
    p = dataCflag(ii); % Index of the node that is cloest to the data point
    x0 = Nd(p,1);
    y0 = Nd(p,2);
    z0 = Nd(p,3);
    if mod(p,length(Y0))~=0 % the node is not on the right edge of the map
       x_hat = Nd(p+1,:) - Nd(p,:);
    else
       x_hat = Nd(p,:) - Nd(p-1,:);
    end
    
    if p > length(Y0) % the node is not on the top edge of the map
        y_hat = Nd(p-length(Y0),:) - Nd(p,:);
    else
        y_hat = Nd(p,:) - Nd(p+length(Y0),:);
    end
    pvector = data(ii,:) - Nd(p,:);
    proj_x = step*( pvector * x_hat' / norm(x_hat) + mod(p,length(Y0))-1 ) + x_lower;
    proj_y = step*( pvector * y_hat' / norm(y_hat) + floor(p/length(Y0)) )+y_lower;
    SOM_proj_data(ii,:) = [proj_x,proj_y];
end

figure(3)
scatter(SOM_proj_data(:,1),SOM_proj_data(:,2));


figure(4)
for i = 1 : N
    if dvar(i) ==1
        plot(SOM_proj_data(i,1), SOM_proj_data(i,2), 'r*');
    elseif dvar(i) ==2
        plot(SOM_proj_data(i,1), SOM_proj_data(i,2), 'ro');
    elseif dvar(i) ==3
        plot(SOM_proj_data(i,1), SOM_proj_data(i,2), 'c+');
    elseif dvar(i) ==4
        plot(SOM_proj_data(i,1), SOM_proj_data(i,2), 'cx');
    elseif dvar(i) ==5
        plot(SOM_proj_data(i,1), SOM_proj_data(i,2), 'bs');
    elseif dvar(i) ==6
        plot(SOM_proj_data(i,1), SOM_proj_data(i,2), 'kd');
    elseif dvar(i) ==7
        plot(SOM_proj_data(i,1), SOM_proj_data(i,2), 'gd');
    end
    hold on;
end



