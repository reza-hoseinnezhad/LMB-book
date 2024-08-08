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