function gz_vals= compute_likelihood(model,z,X)
    M = size(X,1);
    P = X(:,[1 2]);
    Phi= zeros(M,2);
    Phi(:,1) = atan2(P(:,2),P(:,1));
    Phi(:,2) = sqrt(sum(P'.^2))';
    z_error =  repmat(z',[M 1])- Phi;
    D = 3*diag(model.D);
    e_sq = sum( (diag(1./D)*z_error').^2 )';
    gz_vals= exp(-e_sq/2 - log(2*pi*prod(D)));
end
