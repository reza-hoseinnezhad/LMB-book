function Meas = gen_meas(Model, Truth)
    % This function generates measurements based on the Model and Truth.

    % Initialize parameters
    Meas.K = Truth.K;
    Meas.Z = cell(Truth.K, 1);

    % Generate measurements
    for k = 1:Truth.K
        if Truth.N(k) > 0
            % Determine detected target indices
            idx = rand(Truth.N(k), 1) <= compute_pD(Truth.X{k}');
            
            % Generate single target observations if detected
            Meas.Z{k} = gen_observation_fn(Model, Truth.X{k}(:, idx), 'noise');
        end
        
        % Generate clutter points
        N_c = poissrnd(Model.lambda_c); % Number of clutter points
        C = repmat(Model.range_c(:, 1), [1, N_c]) + diag(Model.range_c *...
            [-1; 1]) * rand(Model.z_dim, N_c); % Clutter generation
        
        % Combine detections and clutter
        Meas.Z{k} = [Meas.Z{k} C];
    end
end
