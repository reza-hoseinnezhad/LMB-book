function glmb = lmb2glmb(tt_lmb, H_req)

    % Vector of existence probabilities
    rvect = get_rvals(tt_lmb);
    % Cost vector
    costv = rvect ./ (1 - rvect);
    % Negative log cost
    neglogcostv = -log(costv);
    % K-shortest path to calculate k-best surviving hypotheses/components
    [paths, nlcost] = kshortestwrap_pred(neglogcostv, H_req);

    %paths: It is a cell array with a length equal to the number of required hypotheses/components, i.e., H_req. Each cell in the array contains a vector representing one of the k-shortest paths. The length of each vector in a cell can vary depending on the path. 
    % nlcost: It is a vector with a length equal to the number of required hypotheses/components, i.e., H_req. Each element in the vector represents the negative log cost associated with the corresponding path in the paths cell array.

    % Initialize GLMB distribution
    glmb.tt = tt_lmb;
    
    % Calculate hypotheses/component weights, tracks, and cardinalities
    for hidx = 1:length(nlcost)
        % Weight of hypothesis/component
        glmb.w(hidx) = sum(log(1 - rvect)) - nlcost(hidx); 
        % Tracks in hypothesis/component
        glmb.I{hidx} = paths{hidx}(:);
        % Cardinality of hypothesis/component
        glmb.n(hidx) = length(paths{hidx}(:));
    end
    
    % Normalize weights
    glmb.w = exp(glmb.w - logsumexp(glmb.w));

    % Extract cardinality distribution
    max_card = max(glmb.n);
    for card = 0:max_card
        % Extract probability of n targets
        glmb.cdn(card + 1) = sum(glmb.w(glmb.n == card));
    end
end