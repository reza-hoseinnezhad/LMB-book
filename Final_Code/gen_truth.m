function Truth = gen_truth(Model)
    % This function generates ground truth tracks for multiple targets.
    
    % Initialize parameters
    Truth.K = 100; % length of data/number of scans
    Truth.X = cell(Truth.K, 1); % ground truth for states of targets  
    Truth.N = zeros(Truth.K, 1); % ground truth for number of targets
    Truth.track_list = cell(Truth.K, 1); % absolute index target identities (plotting)
    nbirths = 10; % total number of target births during the simulation
    Truth.total_tracks = nbirths; % total number of appearing tracks
    
    % Define target initial states and birth times
    tbirth = [1, 10, 10, 10, 20, 40, 40, 40, 60, 60];
    tdeath = zeros(nbirths, 1);
    
    xstart = [100, 0, 990, 440, 10, 750, 230, 0, 500, 0;
              10, 300, 550, 1000, 463, 0, 990, 630, 10, 150;
              1, 10, -10, 1, 13, 1, 1, 12, 1, 9;
              10, 1, -1, -10, 1, 13, -10, -1, 9, 1];
    
    % Generate the tracks
    for targetnum = 1:nbirths
        targetstate = xstart(:, targetnum);
        for k = tbirth(targetnum):Truth.K
            Truth.X{k} = [Truth.X{k} targetstate];
            Truth.track_list{k} = [Truth.track_list{k} targetnum];
            Truth.N(k) = Truth.N(k) + 1;
            
            % Update target state
            targetstate = Model.F * targetstate + mvnrnd(zeros(4, 1), Model.Q, 1)';
            
            % Check if the target is out of bounds
            if ~(targetstate(1) >= 0 && targetstate(1) <= 1000 ...
                && targetstate(2) >= 0 && targetstate(2) <= 1000)
                tdeath(targetnum) = k;
                break;
            end
        end
        
        if tdeath(targetnum) == 0
            tdeath(targetnum) = Truth.K;
        end
    end
end
