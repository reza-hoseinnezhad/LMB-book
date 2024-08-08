function [X,W] = gen_newstate_fn(model, Xd, Wd, V)
    % This function generates new states using the given state transition model
    % and an optional noise term.
    %
    % Inputs:
    %   model: model structure containing the state transition matrix (F) and noise covariance (Q)
    %   Xd: input state matrix (dimensions: N x M)
    %   V: noise term (optional); can be either a matrix or the string 'noise' to use model.Q
    %
    % Output:
    %   X: new state matrix (dimensions: M x N)

    % If the input state matrix is empty, return an empty matrix
    if isempty(Xd)
        X = [];
    else
        % Handle the noise term
        if ~isnumeric(V)
            %V = zeros(size(Xd));
            if strcmp(V, 'noise')
                % Generate noise samples from a multivariate normal distribution
                V = mvnrnd(zeros(1,size(Xd,2)), 4*model.Q, size(Xd, 1));
            end
        end

        % Apply the state transition model
        X = (model.F * Xd')' + V;
        idx = (X(:,1)>=0) & (X(:,1)<=1000) ...
                & (X(:,2)>=0) & (X(:,2)<=1000);
        X = X(idx,:);
        W = Wd(idx)*sum(Wd)/sum(Wd(idx)); % picking up only the particle weights for ...
                                            % particles that propagate
                                            % within the area, and
                                            % rescaling them to retain
                                            % their sum so that r
                                            % calculation remains valid.
    end
end
