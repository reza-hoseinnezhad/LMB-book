function pD = compute_pD(X,Sensor_Loc,Model)
    if isempty(X)
        pD= [];
    else
        pD0= 0.98;
        pD_Sigma= 5000;
        M= size(X,1);
        P= X(:,[1 2]) - repmat(Sensor_Loc,[M 1]);
        e_sq= sum((P.^2)')/pD_Sigma^2;
        pD = pD0*exp(-e_sq'/2);
        idx = abs(atan(P(:,2)./P(:,1))) > Model.FOV;
        pD(idx) = eps;
    end
end