function [X,W] = gen_newstate_fn(model, Xd, Wd, V)
    if isempty(Xd)
        X = [];
    else
        % Handle the noise term
        if ~isnumeric(V)
            %V = zeros(size(Xd));
            if strcmp(V, 'noise')
                % Generate noise samples from a 
                % multivariate normal distribution
                V = mvnrnd(zeros(1,size(Xd,2)), 9*model.Q, size(Xd, 1));
            end
        end

        % Apply the state transition model
        X = (model.F * Xd')' + V;
        idx = (X(:,1)>=0) & (X(:,1)<=1000) & (X(:,2)>=0) & (X(:,2)<=1000);
        X = X(idx,:);
        % picking up only the particle weights for particles that propagate within the area, and rescaling them to retain their sum so that r calculation remains valid.
        W = Wd(idx); 
    end
end
