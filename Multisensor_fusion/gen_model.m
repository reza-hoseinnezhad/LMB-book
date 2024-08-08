function model= gen_model
    % basic parameters
    model.x_dim = 4;   %dimension of state vector
    model.z_dim = 2;   %dimension of observation vector
    model.v_dim = 2;   %dimension of process noise
    % dynamical model parameters (CV model)
    % state transformation given by gen_newstate_fn,
    % transition matrix is N/A in non-linear case
    model.T= 1;                        %sampling period
    model.sigma_a = 0.5;
    model.F = [eye(2) model.T*eye(2); zeros(2,2) eye(2)];
    model.Q = model.sigma_a^2*[1/4*model.T^4*eye(2) ... 
        1/2*model.T^3*eye(2); 1/2*model.T^3*eye(2) model.T^2*eye(2)];
    % birth parameters (LMB birth model, single component only)
    model.T_birth = 4;         %no. of LMB birth terms
    %prob of birth for each LMB birth term
    model.r_birth = zeros(model.T_birth,1); 
    model.birthlimits = cell(model.T_birth,1); 
    % First birth term: entry from lower horizontal
    % edge, connecting (0,0)-(1000,0)
    model.r_birth(1) = 0.05;
    model.birthlimits{1} = [0 1000 0 150 -15 15 0.50 15];
    % Second birth term: entry from upper horizontal edge
    model.r_birth(2) = 0.05;
    model.birthlimits{2} = [0 1000 850 1000 -15 15 -0.50 -15];
    % Third birth term: entry from left vertical edge
    model.r_birth(3) = 0.05;
    model.birthlimits{3} = [0 150 0 1000 0.50 15 -15 15];
    % Fourth birth term: entry from right vertical edge
    model.r_birth(4) = 0.05;
    model.birthlimits{4} = [850 1000 0 1000 -15 -0.50 -15 15];
    % std for beaing and range noise detection parameters
    model.D = diag([1*(pi/180); 1]);
    % clutter parameters
    model.lambda_c = 2; %poisson average rate of uniform 
                        % clutter (per scan)
    model.range_c = [0 pi/2; 0 1000*sqrt(2)];          
    model.pdf_c = 1/prod(model.range_c(:,2)-model.range_c(:,1)); 
    %uniform clutter density

    model.sen_loc1 = [0 0];
    model.sen_loc2 = [1000 0];
    model.sen_loc3 = [1000 1000];
    model.sen_loc4 = [0 1000];
    model.FOV = pi/3;
    % setting sensor locations
end