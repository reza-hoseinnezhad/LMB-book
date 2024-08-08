function [tt_lmb_birth, tt_lmb_survive] = lmbpredict(tt_lmb_update, model, filter, k)
    % This function generates birth and surviving tracks for the LMB filter.
    % Inputs:
    %   tt_lmb_update: updated LMB tracks
    %   model: model structure containing parameters and state transition model
    %   filter: filter structure containing filtering parameters
    %   k: current time step

    % Generate birth tracks
    num_birth_tracks = length(model.r_birth);                        % Number of birth tracks
    tt_lmb_birth = cell(num_birth_tracks, 1);                        % Initialize cell array for birth tracks

    for birth_idx = 1:num_birth_tracks
        tt_lmb_birth{birth_idx}.r = model.r_birth(birth_idx);        % Birth probability for birth track
        % Generate samples for birth track
        tt_lmb_birth{birth_idx}.x = gen_4D_rand(model.birthlimits{birth_idx},filter.npt);
        %tt_lmb_birth{birth_idx}.x = mvnrnd(model.m_birth{birth_idx}, model.P_birth{birth_idx}, filter.npt);
        tt_lmb_birth{birth_idx}.w = ones(1,filter.npt)/filter.npt;
        tt_lmb_birth{birth_idx}.l = [k; birth_idx]; % Track label for birth track
    end

    % Generate surviving tracks
    num_surviving_tracks = length(tt_lmb_update);                    % Number of surviving tracks
    tt_lmb_survive = cell(num_surviving_tracks, 1);                  % Initialize cell array for surviving tracks

    for surv_idx = 1:num_surviving_tracks
        % Compute predicted weights for surviving tracks
        wtemp_predict = compute_pS(tt_lmb_update{surv_idx}.x) .* tt_lmb_update{surv_idx}.w(:);

        % Generate predicted samples for surviving tracks
        [xtemp_predict,wtemp_predict] = gen_newstate_fn(model, tt_lmb_update{surv_idx}.x, ...
            wtemp_predict, 'noise');

        % Update surviving track parameters
        tt_lmb_survive{surv_idx}.r = sum(wtemp_predict) * tt_lmb_update{surv_idx}.r;    % Predicted existence probability
        tt_lmb_survive{surv_idx}.x = xtemp_predict;                  % Samples for surviving track
        tt_lmb_survive{surv_idx}.w = wtemp_predict / sum(wtemp_predict);                % Normalized weights
        tt_lmb_survive{surv_idx}.l = tt_lmb_update{surv_idx}.l;      % Track label for surviving track
    end
end
