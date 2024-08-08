function plot_truth(Truth)
    % Define colors for plotting
    num_colors = 10;
    Colors = colormap(lines(num_colors));
    % Create and configure the subplots
    close all;
    set(gcf, 'Position', [3 584 1000 420]);
    subplot(1,2,1);
    set(gca,'XLim',[0 1000],'YLim',[0 1000]);
    hold on;
    subplot(1,2,2);
    set(gca,'XLim',[0 100],'YLim',[0 10]);
    xlabel('time in sec (k)');
    ylabel('cardinality');
    hold on;
    % Plot the ground truth data
    for k = 1:Truth.K
        subplot(1,2,1);
        title(['k = ' num2str(k)]);
        xlabel('x coordinate (m)');
        ylabel('y coordinate (m)');
        if ~isempty(Truth.X{k})
            for i = 1:length(Truth.X{k}(1,:))
                h = plot(Truth.X{k}(1,i),Truth.X{k}(2,i),'.');
                set(h,'Color',Colors(Truth.track_list{k}(i),:));
            end
        end
        subplot(1,2,2);
        plot(k,Truth.N(k),'b.-');
        pause(0.1);
    end
end

