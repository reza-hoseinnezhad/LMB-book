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