function [est,TT_LMB] = run_filter(model,meas)

    est.X= cell(meas.K,1);
    est.N= zeros(meas.K,1);
    est.L= cell(meas.K,1);
    
    %filter parameters
    filter.T_max= 100;                  %maximum number of tracks
    filter.track_threshold= 1e-2;       %threshold to prune tracks
    
    filter.H_bth= 30;                    %requested number of birth components/hypotheses (for LMB to GLMB casting before update)
    filter.H_sur= 1000;                  %requested number of surviving components/hypotheses (for LMB to GLMB casting before update)
    filter.H_upd= 1000;                  %requested number of updated components/hypotheses (for GLMB update)
    filter.H_max= 1000;                  %cap on number of posterior components/hypotheses (not used yet)
    filter.hyp_threshold= 1e-15;         %pruning threshold for components/hypotheses (not used yet)
    
    filter.npt= 1000;                   %number of particles per track
    filter.nth= 100;                    %threshold on effective number of particles before resampling (not used here, resampling is forced at every step, otherwise number of particles per track grows)
    
    filter.run_flag= 'disp';            %'disp' or 'silence' for on the fly output
    
    est.filter= filter;
    
    %=== Filtering
    
    %initial prior
    tt_lmb_update= cell(0,1);      %track table for LMB (cell array of structs for individual tracks)
    TT_LMB = cell(meas.K,1);

    %recursive filtering
    for k=1:meas.K
    
        %prediction
        [tt_lmb_birth,tt_lmb_survive]= lmbpredict(tt_lmb_update,model,filter,k);            
        T_predict= length(tt_lmb_birth)+length(tt_lmb_survive);
    
        %update
        glmb_predict= castlmbpred(tt_lmb_birth,tt_lmb_survive,filter);
        glmb_update= update(glmb_predict,model,filter,meas,k);
        tt_lmb_update= glmb2lmb(glmb_update);                                               
            
        %pruning, truncation and track cleanup
        tt_lmb_update= clean_lmb(tt_lmb_update,filter);                                     
        %T_clean= length(tt_lmb_update);
    
        %state estimation
        [est.X{k},est.N(k),est.L{k}]= extract_estimates(tt_lmb_update,model);
        display(['time = ', num2str(k)]);
        TT_LMB{k,1} = tt_lmb_update;
    end
end
