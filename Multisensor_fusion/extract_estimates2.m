function [X,N,L]=extract_estimates2(tt_lmb,model,est_thresh_r)
    rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999);
    idxcmp = find(rvect>est_thresh_r); 
    N = length(idxcmp);
    X= zeros(model.x_dim,N);
    L= zeros(2,N);
    for n=1:N
        X(:,n)= tt_lmb{idxcmp(n)}.x'*tt_lmb{idxcmp(n)}.w;
        L(:,n)= tt_lmb{idxcmp(n)}.l;
    end
end