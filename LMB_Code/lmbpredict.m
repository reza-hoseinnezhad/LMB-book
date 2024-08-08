function [tt_lmb_birth, tt_lmb_survive] = lmbpredict(tt_lmb_update, model, filter, k)
    % Generate birth tracks

    % Initialize cell array for birth tracks
    tt_lmb_birth = cell(length(model.r_birth), 1);   
                                               
    for birth_idx = 1:length(model.r_birth)
        % Birth probability for birth track
        tt_lmb_birth{birth_idx}.r = model.r_birth(birth_idx); 
        % Generate samples for birth track
        tt_lmb_birth{birth_idx}.x = gen_4D_rand(...
                    model.birthlimits{birth_idx},filter.npt);
        tt_lmb_birth{birth_idx}.w = ones(1,filter.npt)/filter.npt;
        % Track label for birth track
        tt_lmb_birth{birth_idx}.l = [k; birth_idx]; 
    end

    % Generate surviving tracks

    % Number of surviving tracks
    num_surviving_tracks = length(tt_lmb_update);
    % Initialize cell array for surviving tracks
    tt_lmb_survive = cell(num_surviving_tracks, 1);

    for surv_idx = 1:num_surviving_tracks
        % Compute predicted weights for surviving tracks
        wtemp_predict = compute_pS(tt_lmb_update{surv_idx}.x)...
            .* tt_lmb_update{surv_idx}.w(:);

        % Generate predicted samples for surviving tracks
        [xtemp_predict,wtemp_predict] = gen_newstate_fn(model, ...
                tt_lmb_update{surv_idx}.x, wtemp_predict, 'noise');

        % Update surviving track parameters

        % Predicted existence probability
        tt_lmb_survive{surv_idx}.r = sum(wtemp_predict) * ...
            tt_lmb_update{surv_idx}.r;
        % Samples for surviving track
        tt_lmb_survive{surv_idx}.x = xtemp_predict;
        % Normalized weights
        tt_lmb_survive{surv_idx}.w = wtemp_predict / sum(wtemp_predict);
        % Track label for surviving track
        tt_lmb_survive{surv_idx}.l = tt_lmb_update{surv_idx}.l;
    end
end