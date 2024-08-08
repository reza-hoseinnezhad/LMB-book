function glmb_predict = castlmbpred(tt_lmb_birth, tt_lmb_survive, filter)
    % This function generates a predicted GLMB (Generalized Labeled Multi-Bernoulli) 
    % distribution from given birth and surviving LMB (Labeled Multi-Bernoulli) tracks.
    %
    % Inputs:
    %   tt_lmb_birth: cell array of birth LMB tracks
    %   tt_lmb_survive: cell array of surviving LMB tracks
    %   filter: filter structure containing parameters H_bth and H_sur
    %
    % Output:
    %   glmb_predict: predicted GLMB distribution

    % Convert birth and surviving LMBs to GLMB structure
    glmb_birth = lmb2glmb(tt_lmb_birth, filter.H_bth);
    glmb_survive = lmb2glmb(tt_lmb_survive, filter.H_sur);

    % Combine birth and surviving tracks
    glmb_predict.tt = cat(1, glmb_birth.tt, glmb_survive.tt);
    
    % Generate predicted hypotheses/components by convolution of birth and survive GLMBs
    num_birth = length(glmb_birth.w);
    num_survive = length(glmb_survive.w);
    
    for bidx = 1:num_birth
        for sidx = 1:num_survive
            hidx = (bidx - 1) * num_survive + sidx;
            glmb_predict.w(hidx) = glmb_birth.w(bidx) * glmb_survive.w(sidx);
            glmb_predict.I{hidx} = [glmb_birth.I{bidx}; length(glmb_birth.tt) + ...
                                    glmb_survive.I{sidx}];
            glmb_predict.n(hidx) = glmb_birth.n(bidx) + glmb_survive.n(sidx);
        end
    end
    
    % Normalize component weights
    glmb_predict.w = glmb_predict.w / sum(glmb_predict.w);

    % Extract cardinality distribution
    max_card = max(glmb_predict.n);
    for card = 0:max_card
        glmb_predict.cdn(card + 1) = sum(glmb_predict.w(glmb_predict.n == card));
    end
end
