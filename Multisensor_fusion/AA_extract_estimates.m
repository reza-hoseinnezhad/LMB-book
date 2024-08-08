function [X,N,L]=AA_extract_estimates(tt_lmb,model)
    
    N = 0;
    for tabidx = 1:length(tt_lmb)
        if tt_lmb{tabidx}.r > 0.10
            N = N+1;
    end

    X = zeros(model.x_dim,N);
    L = zeros(2,N);
    
    n = 0; 
    for tabidx=1:length(tt_lmb)
        if tt_lmb{tabidx}.r > 0.10
            n = n+1;
            X(:,n)= tt_lmb{tabidx}.x'*tt_lmb{tabidx}.w;
            L(:,n)= tt_lmb{tabidx}.l;
        end
    end
end