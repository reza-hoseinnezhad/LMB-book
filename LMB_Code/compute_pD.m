function pD = compute_pD(X)
    if isempty(X)
        pD= [];
    else
        pD0= 0.98;
        Sensor_Loc= [0 0];
        pD_Sigma= 5000;
        M= size(X,1);
        P= X(:,[1 2]);
        e_sq= sum([(P-repmat(Sensor_Loc,[M 1])).^2]')/pD_Sigma^2;
        pD= pD0*exp(-e_sq'/2);
    end
end