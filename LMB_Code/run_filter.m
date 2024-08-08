function est = run_filter(model,meas)
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
    tt_lmb_update= cell(0,1);   % track table for LMB (cell array of structs for individual tracks)
    %recursive filtering
    for k=1:meas.K
        %prediction
        [tt_lmb_birth,tt_lmb_survive]= lmbpredict(tt_lmb_update, model,filter,k);            
        %update
        glmb_predict= castlmbpred(tt_lmb_birth,tt_lmb_survive,filter);
        glmb_update= update(glmb_predict,model,filter,meas,k);
        tt_lmb_update= glmb2lmb(glmb_update);                                               
        %pruning, truncation and track cleanup
        tt_lmb_update= clean_lmb(tt_lmb_update,filter);                                     
        %state estimation
        [est.X{k},est.N(k),est.L{k}]= extract_estimates(tt_lmb_update,model);
        display(['time = ', num2str(k)]); 
    end
end