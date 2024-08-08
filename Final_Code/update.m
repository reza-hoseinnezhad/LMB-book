function glmb_update = update(glmb_predict, model, filter, meas, k)
    % Update GLMB predictions using measurements
    % Inputs:
    % glmb_predict: prediction struct
    % model: filter model struct
    % filter: filter struct
    % meas: measurement struct
    % k: time index
    % Output:
    % glmb_update: updated GLMB prediction struct
    
    % Get number of measurements
    m = size(meas.Z{k}, 2);
    
    % Initialize cell array for updated tracks
    tt_update = cell((1+m)*length(glmb_predict.tt), 1);
    
    % Missed detection tracks (legacy tracks)
    for tabidx = 1:length(glmb_predict.tt)
        tt_update{tabidx} = glmb_predict.tt{tabidx}; % Same track table
    end
    
    % Measurement updated tracks (all pairs)
    allcostm = zeros(length(glmb_predict.tt), m); % Global cost matrix
    for emm = 1:m
        for tabidx = 1:length(glmb_predict.tt)
            % Index of predicted track i updated with measurement j is
            % (number_predicted_tracks*j + i)
            stoidx = length(glmb_predict.tt)*emm + tabidx;
    
            % Compute weight update for this track and this measurement

            w_temp = compute_pD(glmb_predict.tt{tabidx}.x) .* ...
                glmb_predict.tt{tabidx}.w(:) .* ...
                compute_likelihood(model, meas.Z{k}(:, emm), glmb_predict.tt{tabidx}.x);
    
            % Update particles and weights for the updated track
            x_temp = glmb_predict.tt{tabidx}.x;
            tt_update{stoidx}.x = x_temp;
            tt_update{stoidx}.w = w_temp / sum(w_temp);
            tt_update{stoidx}.l = glmb_predict.tt{tabidx}.l;
    
            % Compute predictive likelihood
            allcostm(tabidx, emm) = sum(w_temp);
        end
    end
    
    % Copy track table back to GLMB struct
    glmb_update.tt = tt_update;
    
    % Precalculation loop for average detection/missed probabilities
    avpd = zeros(length(glmb_predict.tt), 1);
    for tabidx = 1:length(glmb_predict.tt)
        if ~isempty(glmb_predict.tt{tabidx}.w)
            avpd(tabidx) = glmb_predict.tt{tabidx}.w(:)' * compute_pD(glmb_predict.tt{tabidx}.x) + eps(0);
        end
    end
    avqd = 1 - avpd;
    
    % Component updates
    if m == 0
        % No measurements means all missed detections
        for pidx = 1:length(glmb_predict.w)
            glmb_update.w(pidx) = -model.lambda_c + ...
                sum(log(avqd(glmb_predict.I{pidx}))) + ...
                log(glmb_predict.w(pidx)); % Hypothesis/component weight
        end
        glmb_update.I = glmb_predict.I; % Hypothesis/component tracks (via indices to track table)
        glmb_update.n = glmb_predict.n; % Hypothesis/component cardinality
    else
        % Loop over predicted components/hypotheses
        runidx = 1;
        for pidx = 1:length(glmb_predict.w)
            if glmb_predict.n(pidx) == 0
                % No target means all clutter
                glmb_update.w(runidx) = -model.lambda_c + ...
                    m * log(model.lambda_c * model.pdf_c) + ...
                    log(glmb_predict.w(pidx));
                glmb_update.I{runidx} = glmb_predict.I{pidx};
                glmb_update.n(runidx) = glmb_predict.n(pidx);
                runidx = runidx + 1;
            else
                % Perform update for component
                % Calculate best updated hypotheses/components
                costm = allcostm(glmb_predict.I{pidx}, :) ./ ...
                    (model.lambda_c * model.pdf_c * repmat(avqd(glmb_predict.I{pidx}), [1 m])); % Cost matrix
                neglogcostm = -log(costm); % Negative log cost
                [uasses, nlcost] = mbestwrap_updt_custom(neglogcostm, round(filter.H_upd * sqrt(glmb_predict.w(pidx)) / sum(sqrt(glmb_predict.w)))); % Murty's algorithm to calculate m-best assignment hypotheses/components
                % Generate corresponding surviving hypotheses/components
                for hidx = 1:length(nlcost)
                    update_hypcmp_tmp = uasses(hidx, :)';
                    % Hypothesis/component weight
                    glmb_update.w(runidx) = -model.lambda_c + ...
                        m * log(model.lambda_c * model.pdf_c) + ...
                        sum(log(avqd(glmb_predict.I{pidx}))) + ...
                        log(glmb_predict.w(pidx)) - nlcost(hidx);
                    % Hypothesis/component tracks (via indices to track table)
                    glmb_update.I{runidx} = length(glmb_predict.tt) * update_hypcmp_tmp + glmb_predict.I{pidx};
                    % Hypothesis/component cardinality
                    glmb_update.n(runidx) = glmb_predict.n(pidx);
                    runidx = runidx + 1;
                end
            end
        end
    end
    
    % Normalize weights
    glmb_update.w = exp(glmb_update.w - logsumexp(glmb_update.w));
    
    % Extract cardinality distribution
    for card = 0:max(glmb_update.n)
        % Extract probability of n targets
        glmb_update.cdn(card + 1) = sum(glmb_update.w(glmb_update.n == card));
    end
end