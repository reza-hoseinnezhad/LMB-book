function tt_lmb = fused(tt_lmb1,tt_lmb2,tt_lmb3,tt_lmb4)
    % Four LMBs: they have the same number of tracks, with the same labels
    % in the same order. For each track, they have the same particles in
    % the same order. However, some lmbs may have extra particles at the
    % end that can be ignored.
    
    omega1 = 0.25; omega2 = 0.25; omega3 = 0.25; omega4 = 0.25; 
    % assigning equal fusion weights to sensors.

    M = length(tt_lmb1);
    % total number of tracks (same for all)

    %initialize LMB struct
    tt_lmb= cell(M,1);

    for tabidx=1:length(tt_lmb)
       tt_lmb{tabidx}.r = 0;
       tt_lmb{tabidx}.w = [];
       tt_lmb{tabidx}.l = tt_lmb1{tabidx}.l;
       % all sensors return LMB posteriors with the same labeled components
       % so we just take the first set of labels
       tt_lmb{tabidx}.x = tt_lmb1{tabidx}.x;
    end

    for tabidx=1:length(tt_lmb)
         m1 = size(tt_lmb1{tabidx}.x,1);
         m2 = size(tt_lmb2{tabidx}.x,1);
         m3 = size(tt_lmb3{tabidx}.x,1);
         m4 = size(tt_lmb4{tabidx}.x,1);
         m = min([m1 m2 m3 m4]);
         
         r1 = tt_lmb1{tabidx}.r;
         r2 = tt_lmb2{tabidx}.r;
         r3 = tt_lmb3{tabidx}.r;
         r4 = tt_lmb4{tabidx}.r;

         
         w1 = tt_lmb1{tabidx}.w(1:m); 
         if m1>m
             deltam = m1-m;
             if deltam>m, disp([deltam m]); pause; end
             w1(1:deltam) = w1(1:deltam) + tt_lmb1{tabidx}.w((m+1):m1);
             disp(sum(w1));
         end
         
         w2 = tt_lmb2{tabidx}.w(1:m); 
         if m2>m
             deltam = m2-m; 
             if deltam>m, disp([deltam m]); pause; end
             w2(1:deltam) = w2(1:deltam) + tt_lmb2{tabidx}.w((m+1):(m+deltam));
         end

         w3 = tt_lmb3{tabidx}.w(1:m); 
         if m3>m
             deltam = m3-m;
             if deltam>m, disp([deltam m]); pause; end
             w3(1:deltam) = w3(1:deltam) + tt_lmb3{tabidx}.w((m+1):(m+deltam));
         end

         w4 = tt_lmb4{tabidx}.w(1:m); 
         if m4>m
             deltam = m4-m;
             if deltam>m, disp([deltam m]); pause; end
             w4(1:deltam) = w4(1:deltam) + tt_lmb4{tabidx}.w((m+1):(m+deltam));
         end
         
         logw = omega1*log(w1+eps) + omega2*log(w2+eps) +...
                omega3*log(w3+eps) + omega4*log(w4+eps);
         sumw = sum(exp(logw));
         tt_lmb{tabidx}.w = exp(logw)/sumw;
         % fusion of particle weights

         lognum_r = omega1*log(r1+eps) + omega2*log(r2+eps) + ...
                   omega3*log(r3+eps) + omega4*log(r4+eps);
         num_r = exp(lognum_r)*sumw;
         % num_r = sum(((r1*w1).^omega1).*...
         %             ((r2*w2).^omega2).*...
         %             ((r3*w3).^omega3).*...
         %             ((r4*w4).^omega4));

         den_r = num_r + ((1-r1)^omega1)*...
                         ((1-r2)^omega2)*...
                         ((1-r3)^omega3)*...
                         ((1-r4)^omega4);
         tt_lmb{tabidx}.r = num_r/den_r;
         % fusion of existence probabilities
         tt_lmb{tabidx}.x = tt_lmb1{tabidx}.x(1:m,:);
         % replacing the particles with the first m particles of the first
         % lmb, knowing that all the lmbs have the same first m particles.
     end
end