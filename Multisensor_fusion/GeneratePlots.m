function GeneratePlots(Truth,Meas,Estimates)
    close all;

    num_colors = 10;
    Colors = colormap(jet(num_colors));

    K = Meas.K;
    Tracks = cell(Truth.total_tracks,1);
    for k = 1:K
        for i = 1:Truth.total_tracks
            [tf, loc] = ismember(i, Truth.track_list{k});
            if tf
                Tracks{i} = [Tracks{i} Truth.X{k}(:,loc)];
            end
        end
    end
    hold on;
    for i = 1:Truth.total_tracks
        plot(Tracks{i}(1,:),Tracks{i}(2,:),'k');
        h = quiver(Tracks{i}(1,1),Tracks{i}(2,1),10*Tracks{i}(3,1),10*Tracks{i}(4,1),'k','LineWidth',4,'MaxHeadSize',0.8);
    end
    
    % for k = 1:K
    %     Z = Meas.Z{k};
    %     Theta = Z(1,:); Rho = Z(2,:);
    %     plot(Rho.*cos(Theta),Rho.*sin(Theta),'b.');
    % end
    set(gca,'XLim',[0 1000], 'YLim',[0 1000]);
    xlabel('p_x (m)'); ylabel('p_y (m)');

    figure;
    xlabel('p_x (m)'); ylabel('p_y (m)');
    hold on;
    for i = 1:Truth.total_tracks
        plot(Tracks{i}(1,:),Tracks{i}(2,:),'k');
        h = quiver(Tracks{i}(1,1),Tracks{i}(2,1),10*Tracks{i}(3,1),10*Tracks{i}(4,1),'k','LineWidth',4,'MaxHeadSize',0.8);
    end    
    
    Assigned_Labels = []; Assigned_Color = [];
    
    for k = 1:K
        if ~isempty(Estimates.L{k})
            for j = 1:size(Estimates.L{k},2)
                if isempty(Assigned_Labels)
                    i = 1; 
                        Assigned_Color = [Assigned_Color; Colors(i,:)];
                        Assigned_Labels = [Assigned_Labels;Estimates.L{k}(:,j)'];
                else
                    if ~ismember(Estimates.L{k}(:,j)',Assigned_Labels,'rows')
                        i = mod(i,num_colors); i = i+1; 
                        Assigned_Color = [Assigned_Color; Colors(i,:)];
                        Assigned_Labels = [Assigned_Labels; Estimates.L{k}(:,j)'];
                    end
                end
            end
        end
    end

    for k = 1:K
        if ~isempty(Estimates.L{k})
            for j = 1:size(Estimates.L{k},2)
                ell = Estimates.L{k}(:,j)';
                [~,loc] = ismember(ell,Assigned_Labels,'rows');
                DotColor = Assigned_Color(loc,:);
                plot(Estimates.X{k}(1,j),Estimates.X{k}(2,j),'o','MarkerEdgeColor', DotColor,'MarkerFaceColor', DotColor, 'MarkerSize', 4);
            end
        end
    end

    figure;
    plot(Truth.N,'k');
    hold on;
    plot(Estimates.N,'ko');
    xlabel('time (sec)'); ylabel('Cardinality');
end
