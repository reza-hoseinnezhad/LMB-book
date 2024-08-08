function model= gen_model

    % basic parameters
    model.x_dim = 4;   %dimension of state vector
    model.z_dim = 2;   %dimension of observation vector
    model.v_dim = 2;   %dimension of process noise
    
    
    % dynamical model parameters (CV model)
    % state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
    model.T= 1;                        %sampling period
    model.sigma_a = 0.5;
    model.F = [eye(2) model.T*eye(2); zeros(2,2) eye(2)];
    model.Q = model.sigma_a^2*[1/4*model.T^4*eye(2) 1/2*model.T^3*eye(2); 1/2*model.T^3*eye(2) model.T^2*eye(2)];
    
    % survival/death parameters
    model.P_S = 0.99;
    model.Q_S = 1 - model.P_S;
    
    % birth parameters (LMB birth model, single component only)
    model.T_birth = 4;         %no. of LMB birth terms
    
    model.r_birth = zeros(model.T_birth,1); %prob of birth for each LMB birth term
    model.birthlimits = cell(model.T_birth,1); 
    
    % First birth term: entry from lower horizontal edge, connecting (0,0)-(1000,0)
    model.r_birth(1) = 0.05;
    model.birthlimits{1} = [0 1000 0 150 -15 15 0 15];
        
    % Second birth term: entry from upper horizontal edge
    model.r_birth(2) = 0.05;
    model.birthlimits{2} = [0 1000 850 1000 -15 15 -15 0];
    
    
    % Third birth term: entry from left vertical edge
    model.r_birth(3) = 0.05;
    model.birthlimits{3} = [0 150 0 1000 0 15 -15 15];
    
    
    % Fourth birth term: entry from right vertical edge
    model.r_birth(4) = 0.05;
    model.birthlimits{4} = [850 1000 0 1000 -15 0 -15 15];
  
    
    % Fifth (last) birth term: artifical birth for recovering lost tracks
    %model.r_birth(5) = 0.05;
    %model.birthlimits{5} = [50 950 50 950 -50 50 -50 50];
    
    
    % observation model parameters (noisy r/theta only)
    % measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
    model.D = diag([1*(pi/180); 1]);
    model.R = (model.D).^2;      %cov for angle and range noise
    % detection parameters
    model.P_D = .98;   %probability of detection in measurements
    model.Q_D = 1 - model.P_D; %probability of missed detection in measurements
    
    % clutter parameters
    model.lambda_c = 2;                             %poisson average rate of uniform clutter (per scan)
    model.range_c = [0 pi/2; 0 1000*sqrt(2)];          %uniform clutter on r/theta
    model.pdf_c = 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
end