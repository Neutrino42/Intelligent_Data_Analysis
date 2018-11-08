clear all;
close all;

data = load('data.parameters.norm.dat');
var = load('data.dependent_var.norm.dat');
N = size(data,1);
d = size(data,2);
C = cov(data);

% Calculate new base V
[V,D] = eig(C);

eigenvalue = eig(C);

% By default, the "eigenvalue" is in increasing order,
% for visually convience, we want it to be decreasing order.
eigval = zeros(d,1);
for i = 1:d-1
    eigval(i) = eigenvalue(d-i);
end

figure(1);
pareto(eigval);

% Seclect only two base vectors
base = V(:,d-2:d);
base2 = V(:,d-1:d);
% Project the data onto the new axis
Proj_Data = transpose(base) * transpose(data);
Proj_Data_2d = transpose(base2) * transpose(data);



%% Colors string
cl = ['r','g','b','k','y','m','c'];


%% Group the data
x = Proj_Data_2d';
%% Initiate cluster centers (CC)
%set number of clusters
NM = 7;
Qerr=zeros(NM,1);
figure;
for M = 1:NM
     %% Randomise data order
    msize = length(x);
    permidx = randperm(msize);
    x = x(permidx,:);
    for m = 1:M
        CC(m,:) = x(end-m+1,:); % Randomly asign a point x to the codebook vector.
    end
    
    %% First plot
    % Plot the data
    msize = length(x);
    dp = subplot(3,2,1:4);
    plot(x(:,1),x(:,2),'b.','markersize',9); hold on;
    ax = gca;
%    ax.XAxisLocation = 'origin';
 %   ax.YAxisLocation = 'origin';
    axis equal;
    % Plot cluster centers
    plot(CC(:,1),CC(:,2),'*r','markersize',9);
    title('Data Space')
    xp = subplot(3,2,5);
    title('Learning Parameter')
    xlabel ('Iteration Number')
    qp = subplot(3,2,6);
    title('Quantisation Error')
    xlabel ('Number of Clusters')
    %% Start iterations
    it = 1;
    pt = 0.001/M; % pause time
    %% Initialise learning parameter
    tau = 1;
    ksi=ones(6,1)*0;
    ksi(1) = exp(-(it/tau));
    cla(xp)
        subplot(3,2,5);
        bar(ksi,'r');hold on
        plot([0 6],[0.001 0.001],'k--');
        title(['Learning Parameter'])
        xlabel ('Iteration Number')
        set(gca,'YScale','log');
        set(gca,'YLim',[0.0001 1]); hold off
    while ksi(it) > 0.001 % iterations loop
        %% initialise data tag
        dataCflag = zeros(msize,1);
        %% Find distance of all data points to cluster centers
        dists = pdist2(x,CC);
        %% loop through data points
        for ii = 1:msize
            dp = subplot(3,2,1:4);
            plot(x(ii,1),x(ii,2),'dr','markersize',6);
            % Find nearest CC in this point
            [~,p] = min(dists(ii,:));
            plot(CC(p,1),CC(p,2),'or','markersize',9);
            %% Move CC towards this data point
            ve =  x(ii,:)-CC(p,:);%vector between them
            plot ([CC(p,1),x(ii,1)],[CC(p,2),x(ii,2)],'-g','linewidth',2)% plot vector between them
            ve = ve*ksi(it); %vector after ksi aplied
            quiver(CC(p,1),CC(p,2),ve(1),ve(2),0,'-r','MaxHeadSize',9,'linewidth',2) ;hold off;
            CC(p,:) = CC(p,:)+ve; %new CC position
            % update distances for this CC
            cd = pdist2(x,CC(p,:));
            dists(:,p) = cd;
            % update data cluster flag
            dataCflag(ii) = p;
            %% Data Iteration Plot update
            % Plot Data
            pause(pt)
            cla(dp)
            dp = subplot(3,2,1:4);
            plot(x(:,1),x(:,2),'b.','markersize',9); hold on;
            ax = gca;
%            ax.XAxisLocation = 'origin';
 %           ax.YAxisLocation = 'origin';
            axis equal;
            % Plot cluster centers
            plot(CC(:,1),CC(:,2),'*r','markersize',9);
            title('Data Space')
            
             plot(CC(p,1),CC(p,2),'or','markersize',8);
        end
       
        
        title('Data Space')
        %% Plot ksi and quantasation error
        
        
        % Permutate data
        permidx = randperm(msize);
        x = x(permidx,:);
        pt = pt*0.1; % increase plot speed for each iteration
        it = it + 1; % increase iteration
        %% Update Learning parameter
        ksi(it) = exp(-(it/tau));
        cla(xp)
        subplot(3,2,5);
        bar(ksi,'r');hold on
        plot([0 7.5],[0.001 0.001],'k--');
        title(['Learning Parameter'])
        xlabel ('Iteration Number')
        set(gca,'YScale','log');
        set(gca,'YLim',[0.0001 1]); hold off
    end
    %% Plot clusters for this number of clusters
        cla(dp)
        dp = subplot(3,2,1:4);
        dataCflag = dataCflag(permidx);
        for cn = 1:M
            plot(x(dataCflag==cn,1),x(dataCflag==cn,2),[cl(cn),'.']); hold on;
            plot(CC(cn,1),CC(cn,2),[cl(cn),'o'],'markersize',8);
            plot(CC(cn,1),CC(cn,2),[cl(cn),'x'],'markersize',8);
            ax = gca;
   %         ax.XAxisLocation = 'origin';
   %         ax.YAxisLocation = 'origin';
            axis equal;
        end
        hold off;
     %% Calculate quantasation error for this number of clusters
        for cn = 1:M
            Qerr(M) = Qerr(M)+ sum(pdist2(x(dataCflag==cn,:),CC(cn,:)));
        end
        Qerr(M)=Qerr(M)/msize;
        cla(qp)
        qp = subplot(3,2,6);
        bar(Qerr);
        set(gca,'YLim',[0 Qerr(1)]);
        xlabel ('Number of Clusters')
        title(['Quantisation Error'])
        pause(4)
        cla(xp)
        subplot(3,2,5);
        title(['Learning Parameter'])
end





figure(3);
for i = 1 : N
    if var(i) ==1
        plot(Proj_Data_2d(1,i), Proj_Data_2d(2,i), 'r*');
    elseif var(i) ==2
        plot(Proj_Data_2d(1,i), Proj_Data_2d(2,i), 'ro');
    elseif var(i) ==3
        plot(Proj_Data_2d(1,i), Proj_Data_2d(2,i), 'c+');
    elseif var(i) ==4
        plot(Proj_Data_2d(1,i), Proj_Data_2d(2,i), 'cx');
    elseif var(i) ==5
        plot(Proj_Data_2d(1,i), Proj_Data_2d(2,i), 'bs');
    elseif var(i) ==6
        plot(Proj_Data_2d(1,i), Proj_Data_2d(2,i), 'kd');
    elseif var(i) ==7
        plot(Proj_Data_2d(1,i), Proj_Data_2d(2,i), 'gd');
    end
    hold on;
end
hold off;


figure(4);

for i = 1 : N
    if var(i) ==1
        plot3(Proj_Data(1,i), Proj_Data(2,i), Proj_Data(3,i), 'r*');
    elseif var(i) ==2
        plot3(Proj_Data(1,i), Proj_Data(2,i), Proj_Data(3,i), 'ro');
    elseif var(i) ==3
        plot3(Proj_Data(1,i), Proj_Data(2,i), Proj_Data(3,i), 'c+');
    elseif var(i) ==4
        plot3(Proj_Data(1,i), Proj_Data(2,i), Proj_Data(3,i), 'cx');
    elseif var(i) ==5
        plot3(Proj_Data(1,i), Proj_Data(2,i), Proj_Data(3,i), 'bs');
    elseif var(i) ==6
        plot3(Proj_Data(1,i), Proj_Data(2,i), Proj_Data(3,i), 'kd');
    elseif var(i) ==7
        plot3(Proj_Data(1,i), Proj_Data(2,i), Proj_Data(3,i), 'gd');
    end
    hold on;
end
hold off;



