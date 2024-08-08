function [X,N,L]=extract_estimates(tt_lmb,model)
    rvect= get_rvals(tt_lmb); rvect= min(rvect,0.999);
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