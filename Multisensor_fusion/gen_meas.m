function Meas = gen_meas(Model, Truth)
    % Initialize parameters
    Meas.K = Truth.K;
    Meas.Z1 = cell(Truth.K, 1);
    Meas.Z2 = cell(Truth.K, 1);
    Meas.Z3 = cell(Truth.K, 1);
    Meas.Z4 = cell(Truth.K, 1);
    % Generate measurements
    for k = 1:Truth.K
        if Truth.N(k) > 0
            % Determine detected target indices
            idx = rand(Truth.N(k), 1) <= compute_pD(Truth.X{k}',[0 0],Model);
            % Generate single target observations if detected
            Meas.Z1{k} = gen_observation_fn(Model, Truth.X{k}(:, idx),Model.sen_loc1,'noise');
            idx = rand(Truth.N(k), 1) <= compute_pD(Truth.X{k}',[1000 0],Model);
            Meas.Z2{k} = gen_observation_fn(Model, Truth.X{k}(:, idx),Model.sen_loc2,'noise');
            idx = rand(Truth.N(k), 1) <= compute_pD(Truth.X{k}',[1000 1000],Model);
            Meas.Z3{k} = gen_observation_fn(Model, Truth.X{k}(:, idx),Model.sen_loc3,'noise');
            idx = rand(Truth.N(k), 1) <= compute_pD(Truth.X{k}',[0 1000],Model);
            Meas.Z4{k} = gen_observation_fn(Model, Truth.X{k}(:, idx),Model.sen_loc4,'noise');
        end
        % Generate clutter points
        N_c1 = poissrnd(Model.lambda_c); % Number of clutter points
        C1 = repmat(Model.range_c(:, 1), [1, N_c1]) + diag(Model.range_c *[-1; 1]) * rand(Model.z_dim, N_c1);
        idx = abs(atan(C1(2,:)./C1(1,:)))< Model.FOV;
        C1 = C1(:,idx);
        N_c2 = poissrnd(Model.lambda_c); % Number of clutter points
        C2 = repmat(Model.range_c(:, 1), [1, N_c2]) + diag(Model.range_c *[-1; 1]) * rand(Model.z_dim, N_c2);
        idx = abs(atan(C2(2,:)./(1000-C2(1,:))))< Model.FOV;
        C2 = C2(:,idx);
        N_c3 = poissrnd(Model.lambda_c); % Number of clutter points
        C3 = repmat(Model.range_c(:, 1), [1, N_c3]) + diag(Model.range_c *[-1; 1]) * rand(Model.z_dim, N_c3);
        idx = abs(atan((1000-C3(2,:))./(1000-C3(1,:))))< Model.FOV;
        C3 = C3(:,idx);
        N_c4 = poissrnd(Model.lambda_c); % Number of clutter points
        C4 = repmat(Model.range_c(:, 1), [1, N_c4]) + diag(Model.range_c *[-1; 1]) * rand(Model.z_dim, N_c4);
        idx = abs(atan((1000-C4(2,:))./C4(1,:)))< Model.FOV;
        C4 = C4(:,idx);
        

        % Combine detections and clutter
        Meas.Z1{k} = [Meas.Z1{k} C1];
        Meas.Z2{k} = [Meas.Z2{k} C2];
        Meas.Z3{k} = [Meas.Z3{k} C3];
        Meas.Z4{k} = [Meas.Z4{k} C4];
    end
end

