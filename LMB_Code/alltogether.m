Model = gen_model;
Truth = gen_truth(Model);
plot_truth(Truth);
Meas=  gen_meas(Model,Truth);
[Estimates,LMBs] =   run_filter(Model,Meas);
GeneratePlots(Truth,Meas,Estimates);

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
    % survival/death parameters
    model.P_S = 0.99;
    model.Q_S = 1 - model.P_S;
    % birth parameters (LMB birth model, single component only)
    model.T_birth = 4;         %no. of LMB birth terms
    model.r_birth = zeros(model.T_birth,1); %prob of birth for 
    % each LMB birth term
    model.birthlimits = cell(model.T_birth,1); 
    % First birth term: entry from lower horizontal ...
    % edge, connecting (0,0)-(1000,0)
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
    % observation model parameters (noisy r/theta only)
    % measurement transformation given by gen_observation_fn,
    % observation matrix is N/A in non-linear case
    model.D = diag([1*(pi/180); 1]); %std for angle and range noise 
    % detection parameters
    model.P_D = .98;   %probability of detection in measurements
    model.Q_D = 1 - model.P_D; %probability of missed detection 
    % in measurements
    % clutter parameters
    model.lambda_c = 2; %poisson average rate of uniform 
    % clutter (per scan)
    model.range_c = [0 pi/2; 0 1000*sqrt(2)];          
    model.pdf_c = 1/prod(model.range_c(:,2)-model.range_c(:,1)); 
    %uniform clutter density
end

function plot_truth(Truth)
    % Define colors for plotting
    num_colors = 10;
    Colors = colormap(lines(num_colors));
    % Create and configure the subplots
    close all;
    set(gcf, 'Position', [3 584 1000 420]);
    subplot(1,2,1);
    set(gca,'XLim',[0 1000],'YLim',[0 1000]);
    hold on;track_list
    subplot(1,2,2);
    set(gca,'XLim',[0 100],'YLim',[0 10]);
    xlabel('time in sec (k)');
    ylabel('cardinality');
    hold on;
    % Plot the ground truth data
    for k = 1:Truth.K
        subplot(1,2,1);
        title(['k = ' num2str(k)]);
        xlabel('x coordinate (m)');
        ylabel('y coordinate (m)');
        if ~isempty(Truth.X{k})
            for i = 1:length(Truth.X{k}(1,:))
                h = plot(Truth.X{k}(1,i),Truth.X{k}(2,i),'.');
                set(h,'Color',Colors(Truth.track_list{k}(i),:));
            end
        end
        subplot(1,2,2);
        plot(k,Truth.N(k),'b.-');
        pause(0.1);
    end
end

function Meas = gen_meas(Model, Truth)
    % Initialize parameters
    Meas.K = Truth.K;
    Meas.Z = cell(Truth.K, 1);
    % Generate measurements
    for k = 1:Truth.K
        if Truth.N(k) > 0
            % Determine detected target indices
            idx = rand(Truth.N(k), 1) <= compute_pD(Truth.X{k}');
            % Generate single target observations if detected
            Meas.Z{k} = gen_observation_fn(Model, ...
                Truth.X{k}(:, idx), 'noise');
        end
        % Generate clutter points
        N_c = poissrnd(Model.lambda_c); % Number ...
        % of clutter points
        C = repmat(Model.range_c(:, 1), [1, N_c]) ...
            + diag(Model.range_c *[-1; 1]) * ...
            rand(Model.z_dim, N_c); % Clutter generation
        % Combine detections and clutter
        Meas.Z{k} = [Meas.Z{k} C];
    end
end


function pD = compute_pD(X)
    if isempty(X)
        pD= [];
    else
        pD0= 0.98;
        Sensor_Loc= [0 0];
        pD_Sigma= 5000;
        M= size(X,1);
        P= X(:,[1 2]);
        e_sq= sum([(P-repmat(Sensor_Loc,[M 1])).^2]')/pD_Sigma^2;
        pD= pD0*exp(-e_sq'/2);
    end
end

function Z= gen_observation_fn(model,X,W)
    if ~isnumeric(W)
        if strcmp(W,'noise')
            W= model.D*randn(size(model.D,2),size(X,2));
        elseif strcmp(W,'noiseless')
            W= zeros(size(model.D,1),size(X,2));
        end
    end
    if isempty(X)
        Z= [];
    else 
        P= X([1 2],:);
        Z(1,:)= atan2(P(2,:),P(1,:));   
        Z(2,:)= sqrt(sum(P.^2));
        Z= Z+ W;
    end
end

function [est,TT_LMB] = run_filter(model,meas)
    est.X= cell(meas.K,1);
    est.N= zeros(meas.K,1);
    est.L= cell(meas.K,1);
    %filter parameters
    filter.T_max= 100;                  %maximum number of tracks
    filter.track_threshold= 1e-2;       %threshold to prune tracks
    filter.H_bth= 30;           %requested number of birth components
                                % or hypotheses (for LMB to GLMB 
                                % casting before update)
    filter.H_sur= 1000;         % requested number of surviving components
                                % or hypotheses (for LMB to GLMB casting 
                                % before update)
    filter.H_upd= 1000;         % requested number of updated components
                                % or hypotheses (for GLMB update)
    filter.H_max= 1000;         % cap on number of posterior components
                                % or hypotheses (not used yet)
    filter.npt= 1000;           % number of particles per track
    est.filter= filter;
    
    %=== Filtering
    
    %initial prior
    tt_lmb_update= cell(0,1);   % track table for LMB (cell array 
                                % of structs for individual tracks)
    TT_LMB = cell(meas.K,1);
    %recursive filtering
    for k=1:meas.K
        %prediction
        [tt_lmb_birth,tt_lmb_survive]= lmbpredict(tt_lmb_update,...
            model,filter,k);            
        T_predict= length(tt_lmb_birth)+length(tt_lmb_survive);
        %update
        glmb_predict= castlmbpred(tt_lmb_birth,tt_lmb_survive,filter);
        glmb_update= update(glmb_predict,model,filter,meas,k);
        tt_lmb_update= glmb2lmb(glmb_update);                                               
        %pruning, truncation and track cleanup
        tt_lmb_update= clean_lmb(tt_lmb_update,filter);                                     
        %state estimation
        [est.X{k},est.N(k),est.L{k}]= extract_estimates(tt_lmb_update...
                                        ,model);
        display(['time = ', num2str(k)]);
        TT_LMB{k,1} = tt_lmb_update;
    end
end

function [tt_lmb_birth, tt_lmb_survive] = lmbpredict(tt_lmb_update, ...
                                            model, filter, k)
    % Generate birth tracks
    num_birth_tracks = length(model.r_birth);   % Number of birth tracks
    tt_lmb_birth = cell(num_birth_tracks, 1);   % Initialize cell array 
                                                % for birth tracks

    for birth_idx = 1:num_birth_tracks
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

function Particles = gen_4D_rand(limits,npt)
    Particles = zeros(npt,4);
    xmin = limits(1);     xmax = limits(2);
    ymin = limits(3);     ymax = limits(4);
    xdotmin = limits(5);     xdotmax = limits(6);
    ydotmin = limits(7);     ydotmax = limits(8);
    Particles(:,1) = (xmax - xmin)*rand(npt,1) + xmin;
    Particles(:,2) = (ymax - ymin)*rand(npt,1) + ymin;
    Particles(:,3) = (xdotmax - xdotmin)*rand(npt,1) + xdotmin;
    Particles(:,4) = (ydotmax - ydotmin)*rand(npt,1) + ydotmin;
end

function pS = compute_pS(X)
    if isempty(X)
        pS= [];
    else
        pS= 0.99*ones(size(X,1),1);
        pS= pS(:);
    end
end

function [X,W] = gen_newstate_fn(model, Xd, Wd, V)
    if isempty(Xd)
        X = [];
    else
        % Handle the noise term
        if ~isnumeric(V)
            %V = zeros(size(Xd));
            if strcmp(V, 'noise')
                % Generate noise samples from a 
                % multivariate normal distribution
                V = mvnrnd(zeros(1,size(Xd,2)), 4*model.Q, size(Xd, 1));
            end
        end

        % Apply the state transition model
        X = (model.F * Xd')' + V;
        idx = (X(:,1)>=0) & (X(:,1)<=1000) ...
                & (X(:,2)>=0) & (X(:,2)<=1000);
        X = X(idx,:);
        % picking up only the particle weights for particles 
        % that propagate within the area, and rescaling them to 
        % retain their sum so that r calculation remains valid.
        W = Wd(idx)*sum(Wd)/sum(Wd(idx)); 
    end
end

function glmb_predict = castlmbpred(tt_lmb_birth, tt_lmb_survive, filter)

    % Convert birth and surviving LMBs to GLMB structure
    glmb_birth = lmb2glmb(tt_lmb_birth, filter.H_bth);
    glmb_survive = lmb2glmb(tt_lmb_survive, filter.H_sur);

    % Combine birth and surviving tracks
    glmb_predict.tt = cat(1, glmb_birth.tt, glmb_survive.tt);
    
    % Generate predicted hypotheses/components by convolution
    % of birth and survive GLMBs
    num_birth = length(glmb_birth.w);
    num_survive = length(glmb_survive.w);
    
    for bidx = 1:num_birth
        for sidx = 1:num_survive
            hidx = (bidx - 1) * num_survive + sidx;
            glmb_predict.w(hidx) = glmb_birth.w(bidx) * ...
                                    glmb_survive.w(sidx);
            glmb_predict.I{hidx} = [glmb_birth.I{bidx}; ...
                    length(glmb_birth.tt) + glmb_survive.I{sidx}];
            glmb_predict.n(hidx) = glmb_birth.n(bidx) + ...
                                    glmb_survive.n(sidx);
        end
    end
    
    % Normalize component weights
    glmb_predict.w = glmb_predict.w / sum(glmb_predict.w);

    % Extract cardinality distribution
    max_card = max(glmb_predict.n);
    for card = 0:max_card
        glmb_predict.cdn(card + 1) = ...
            sum(glmb_predict.w(glmb_predict.n == card));
    end
end

function glmb = lmb2glmb(tt_lmb, H_req)

    % Vector of existence probabilities
    rvect = get_rvals(tt_lmb);
    % Cost vector
    costv = rvect ./ (1 - rvect);
    % Negative log cost
    neglogcostv = -log(costv);
    % K-shortest path to calculate k-best surviving hypotheses/components
    [paths, nlcost] = kshortestwrap_pred(neglogcostv, H_req);

    %paths: It is a cell array with a length equal to the number of 
    % required hypotheses/components, i.e., H_req. Each cell in the 
    % array contains a vector representing one of the k-shortest paths
    % The length of each vector in a cell can vary depending on the path.
    % 
    % nlcost: It is a vector with a length equal to the number of
    % required hypotheses/components, i.e., H_req. Each element in
    % the vector represents the negative log cost associated with
    % the corresponding path in the paths cell array.

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

function rvect= get_rvals(tt_lmb)
% function to extract vector of existence probabilities 
% from LMB track table
    rvect= zeros(length(tt_lmb),1);
    for tabidx=1:length(tt_lmb)
        rvect(tabidx)= tt_lmb{tabidx}.r;
    end
end

function logsum = logsumexp(w)
    %performs log-sum-exp trick to avoid numerical underflow
    %input:  w weight vector assumed already log transformed
    %output: log(sum(exp(w)))
    if all(w==-inf)
        logsum= -inf;
        return;
    end
    
    [val,~] = max(w);
    logsum = log(sum(exp(w-val))) + val;
end

function glmb_update = update(glmb_predict, model, filter, meas, k)
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
                compute_likelihood(model, meas.Z{k}(:, emm), ...
                glmb_predict.tt{tabidx}.x);
    
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
            avpd(tabidx) = glmb_predict.tt{tabidx}.w(:)' ...
                * compute_pD(glmb_predict.tt{tabidx}.x) + eps(0);
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
        glmb_update.I = glmb_predict.I; % Hypothesis/component tracks
                                        % (via indices to track table)
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
                    (model.lambda_c * model.pdf_c * ...
                    repmat(avqd(glmb_predict.I{pidx}), ...
                    [1 m])); % Cost matrix
                neglogcostm = -log(costm); % Negative log cost
                % Murty's algorithm to calculate m-best 
                % assignment hypotheses/components
                [uasses, nlcost] = mbestwrap_updt_custom(...
                    neglogcostm, round(filter.H_upd * ...
                    sqrt(glmb_predict.w(pidx)) / ...
                    sum(sqrt(glmb_predict.w)))); 
                % Generate corresponding surviving hypotheses
                % or components
                for hidx = 1:length(nlcost)
                    update_hypcmp_tmp = uasses(hidx, :)';
                    % Hypothesis/component weight
                    glmb_update.w(runidx) = -model.lambda_c + ...
                        m * log(model.lambda_c * model.pdf_c) + ...
                        sum(log(avqd(glmb_predict.I{pidx}))) + ...
                        log(glmb_predict.w(pidx)) - nlcost(hidx);
                    % Hypothesis/component tracks (via 
                    % indices to track table)
                    glmb_update.I{runidx} = length(glmb_predict.tt)...
                        * update_hypcmp_tmp + glmb_predict.I{pidx};
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
        glmb_update.cdn(card + 1) = ...
            sum(glmb_update.w(glmb_update.n == card));
    end
end

function gz_vals= compute_likelihood(model,z,X)
    M= size(X,1);
    P= X(:,[1 2]);
    Phi= zeros(M,2);
    Phi(:,1)= atan2(P(:,2),P(:,1));
    Phi(:,2)= sqrt(sum(P'.^2))';
    z_error =  repmat(z',[M 1])- Phi;
    D = diag(model.D);
    e_sq= sum( (diag(1./D)*z_error').^2 )';
    gz_vals= exp(-e_sq/2 - log(2*pi*prod(D)));
end

function tt_lmb= glmb2lmb(glmb)

    % find unique labels (with different possibly 
    % different association histories)
    lmat= zeros(2,length(glmb.tt),1); 
    for tabidx= 1:length(glmb.tt)
        lmat(:,tabidx)= glmb.tt{tabidx}.l;
    end
    lmat= lmat'; %matrix of labels, each label is a row
    %sort(lmat);
    
    [cu,~,ic]= unique(lmat,'rows'); cu= cu'; %cu is a matrix,
                                             % each label is a column
    
    %initialize LMB struct
    tt_lmb= cell(size(cu,2),1);
    
    for tabidx=1:length(tt_lmb)
       tt_lmb{tabidx}.r= 0;
       tt_lmb{tabidx}.x= [];
       tt_lmb{tabidx}.w= []; %particle weights
       tt_lmb{tabidx}.l= cu(:,tabidx);
    end
    
    Num_particles = zeros(1,length(tt_lmb));
    Particle_Index = ones(1,length(tt_lmb));
    
    for hidx=1:length(glmb.w) %for each hypothesis  
       for t= 1:glmb.n(hidx)    %for each element in the hypothesis
          trkidx= glmb.I{hidx}(t);  %find the t'th track
          newidx= ic(trkidx);
          Num_particles(newidx)  = Num_particles(newidx) + ...
              length(glmb.tt{trkidx}.w); %num particles in
                                         % each new lmb component
       end
    end
    
    for tabidx=1:length(tt_lmb)
       tt_lmb{tabidx}.x= zeros(Num_particles(tabidx),4);
       tt_lmb{tabidx}.w= zeros(Num_particles(tabidx),1);
    end
    
    
    %extract individual tracks
    for hidx=1:length(glmb.w) %each hypothesis
       for t= 1:glmb.n(hidx) %each element in a hypothesis
          trkidx= glmb.I{hidx}(t); %element t of hypothesis 
          newidx= ic(trkidx);
                
          tt_lmb{newidx}.x(Particle_Index(newidx):...
              (Particle_Index(newidx)+length(glmb.tt{trkidx}.w)-1),:)...
              = glmb.tt{trkidx}.x;
          tt_lmb{newidx}.w(Particle_Index(newidx):...
              (Particle_Index(newidx)+length(glmb.tt{trkidx}.w)-1)) ...
              = glmb.w(hidx)*glmb.tt{trkidx}.w;
          Particle_Index(newidx) = Particle_Index(newidx) + ...
              length(glmb.tt{trkidx}.w);
       end
    end
    
    %extract existence probabilities and normalize track weights
    for tabidx=1:length(tt_lmb)
       tt_lmb{tabidx}.r= sum(tt_lmb{tabidx}.w); % ensures no ...
       % potential matlab rounding error to 1
       tt_lmb{tabidx}.w= tt_lmb{tabidx}.w/tt_lmb{tabidx}.r;
    end
end

function tt_lmb_out= clean_lmb(tt_lmb_in,filter)
    %prune tracks with low existence probabilities

    rvect= get_rvals(tt_lmb_in);
    idxkeep= find(rvect > filter.track_threshold);
    rvect= rvect(idxkeep);
    tt_lmb_out= tt_lmb_in(idxkeep);
    
    %enforce cap on maximum number of tracks
    if length(tt_lmb_out) > filter.T_max
        [~,idxkeep]= sort(rvect,'descend');
        tt_lmb_out= tt_lmb_out(idxkeep(1:filter.T_max));
    end
    
    %resampling the particles in each tracks
    for tabidx=1:length(tt_lmb_out)
        xtemptemp= tt_lmb_out{tabidx}.x;
        wtemptemp= tt_lmb_out{tabidx}.w;
        rspidx= randsample(length(wtemptemp),filter.npt,true,wtemptemp); 
        tt_lmb_out{tabidx}.x= xtemptemp(rspidx,:);
        tt_lmb_out{tabidx}.w= ones(filter.npt,1)/filter.npt;
    end
end

function [X,N,L]=extract_estimates(tt_lmb,model)
    rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999); 
    rvect= max(rvect,0.001);
    cdn= prod(1-rvect)*esf(rvect./(1-rvect));
    [~,mode] = max(cdn);
    N = min(length(rvect),mode-1);
    X= zeros(model.x_dim,N);
    L= zeros(2,N);
    
    [~,idxcmp]= sort(rvect,'descend');
    for n=1:N
        X(:,n)= tt_lmb{idxcmp(n)}.x'*tt_lmb{idxcmp(n)}.w;
        L(:,n)= tt_lmb{idxcmp(n)}.l;
    end
end