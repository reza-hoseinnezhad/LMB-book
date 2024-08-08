function Z= gen_observation_fn(model,X,sen_loc,W)
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
        P= X([1 2],:)-sen_loc';
        Z(1,:)= abs(atan(P(2,:)./P(1,:)));   
        Z(2,:)= sqrt(sum(P.^2));
        idx = ( Z(1,:)<model.FOV );
        Z= Z(:,idx)+ W(:,idx);
    end
end