function pS = compute_pS(X)

if isempty(X)
    pS= [];
else
    pS= 0.99*ones(size(X,1),1);
    pS= pS(:);
end
