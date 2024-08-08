function gz_vals= compute_likelihood(model,z,X)

% compute likelihood vector g= [ log_g(z|x_1), ... , log_g(z|x_M) ] -
% this is for bearings and range case with additive Gaussian noise
M= size(X,1);
P= X(:,[1 2]);
Phi= zeros(M,2);
Phi(:,1)= atan2(P(:,1),P(:,2));
Phi(:,2)= sqrt(sum(P'.^2))';
z_error =  repmat(z',[M 1])- Phi;
D = 10*diag(model.D);
e_sq= sum( (diag(1./D)*z_error').^2 )';
gz_vals= exp(-e_sq/2 - log(2*pi*prod(D)));
