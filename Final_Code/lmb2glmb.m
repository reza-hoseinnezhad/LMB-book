function glmb = lmb2glmb(tt_lmb, H_req)
    % This function converts an LMB (Labeled Multi-Bernoulli) distribution
    % to a GLMB (Generalized Labeled Multi-Bernoulli) distribution.
    %
    % Inputs:
    %   tt_lmb: cell array of LMB tracks
    %   H_req: number of required hypotheses/components
    %
    % Output:
    %   glmb: GLMB distribution

    % Vector of existence probabilities
    rvect = get_rvals(tt_lmb);
    
    % Cost vector
    costv = rvect ./ (1 - rvect);
    
    % Negative log cost
    neglogcostv = -log(costv);
    
    % K-shortest path to calculate k-best surviving hypotheses/components
    [paths, nlcost] = kshortestwrap_pred(neglogcostv, H_req);

    %paths: It is a cell array with a length equal to the number of ...
    % required hypotheses/components, i.e., H_req. Each cell in the ...
    % array contains a vector representing one of the k-shortest paths...
    % The length of each vector in a cell can vary depending on the path.
    % 
    % nlcost: It is a vector with a length equal to the number of ...
    % required hypotheses/components, i.e., H_req. Each element in ...
    % the vector represents the negative log cost associated with ... 
    % the corresponding path in the paths cell array.

    % Initialize GLMB distribution
    glmb.tt = tt_lmb;
    
    % Calculate hypotheses/component weights, tracks, and cardinalities
    for hidx = 1:length(nlcost)
        glmb.w(hidx) = sum(log(1 - rvect)) - nlcost(hidx); % Weight of hypothesis/component
        glmb.I{hidx} = paths{hidx}(:);                      % Tracks in hypothesis/component
        glmb.n(hidx) = length(paths{hidx}(:));              % Cardinality of hypothesis/component
    end
    
    % Normalize weights
    glmb.w = exp(glmb.w - logsumexp(glmb.w));

    % Extract cardinality distribution
    max_card = max(glmb.n);
    for card = 0:max_card
        glmb.cdn(card + 1) = sum(glmb.w(glmb.n == card)); % Extract probability of n targets
    end
end