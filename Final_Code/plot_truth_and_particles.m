function plot_truth_and_particles(Truth,LMB)
    figure;
    set(gca,'XLim',[0 1000],'YLim',[0 1000]);
    xlabel('x (m)');
    ylabel('y(m)');
    num_colors = 10;
    Colors = colormap(jet(num_colors));

    for k = 1:Truth.K
        title([' k = ', num2str(k)]);
        hold on;
        True_X = Truth.X{k}; 
        for i = 1:size(True_X,2)
            plot(True_X(1,i),True_X(2,i),'rO');
        end

        if ~isempty(LMB{k,1})
            for i = 1:length(LMB{k,1})
                xtemp = LMB{k,1}(i);
                particles = xtemp{1}.x;
                plot(particles(:,1),particles(:,2),'.','Color',Colors(mod(i,10)+1,:));
                xhat = mean(particles(:,1));
                yhat = mean(particles(:,2));
                text(xhat,yhat,['r = ',num2str(xtemp{1}.r)]);
            end
        end
        pause;
        hold off;
        delete(gca);
        set(gca,'XLim',[0 1000],'YLim',[0 1000]);
    end
end