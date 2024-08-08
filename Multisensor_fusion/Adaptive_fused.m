function tt_lmb = Adaptive_fused(tt_lmb1,tt_lmb2,tt_lmb3,tt_lmb4)
    % Four LMBs: they have the same number of tracks, with the same labels in the same order. For each track, they have the same particles in the same order. However, there are reptitions of particles with the particle weights being spread through the repetitions.
    
    M = length(tt_lmb1);
    % total number of tracks (same for all)
    tt_lmb= cell(M,1);
    %initialize LMB struct
    for tabidx=1:length(tt_lmb)
       tt_lmb{tabidx}.r = 0; tt_lmb{tabidx}.w = []; tt_lmb{tabidx}.l = tt_lmb1{tabidx}.l;
       % all sensors return LMB posteriors with the same labeled components so we just take the first set of labels
    end
    for tabidx=1:length(tt_lmb)
         r1 = tt_lmb1{tabidx}.r;
         r2 = tt_lmb2{tabidx}.r;
         r3 = tt_lmb3{tabidx}.r;
         r4 = tt_lmb4{tabidx}.r;
         [tt_lmb{tabidx}.x, ~, ic] = unique(tt_lmb1{tabidx}.x,'rows','stable');  
         % Find unique elements in tt_lmb1{tabidx}.x, keep order
         w1 = accumarray(ic, tt_lmb1{tabidx}.w); w1 = w1/sum(w1); 
         % Sum up elements in tt_lmb1{tabidx}.w corresponding to each unique 
         % element in tt_lmb{tabidx}.x
         [~, ~, ic] = unique(tt_lmb2{tabidx}.x,'rows','stable');
         w2 = accumarray(ic, tt_lmb2{tabidx}.w); w2 = w2/sum(w2);
         [~, ~, ic] = unique(tt_lmb3{tabidx}.x,'rows','stable');
         w3 = accumarray(ic, tt_lmb3{tabidx}.w); w3 = w3/sum(w3);
         [~, ~, ic] = unique(tt_lmb4{tabidx}.x,'rows','stable');
         w4 = accumarray(ic, tt_lmb4{tabidx}.w); w4 = w4/sum(w4);

         % assigning adaptive fusion weights to sensors.
         sumr = r1 + r2 + r3 + r4;
         omega1 = r1/sumr; omega2 = r2/sumr; omega3 = r3/sumr; omega4 = r4/sumr; 
         
         rw = ((r1*w1).^omega1) .* ((r2*w2).^omega2) .* ((r3*w3).^omega3) .* ((r4*w4).^omega4); 
         tt_lmb{tabidx}.w = rw/sum(rw);
         % fusion of particle weights
         den_r = ((1-r1)^omega1)*((1-r2)^omega2)*((1-r3)^omega3)*((1-r4)^omega4) + sum(rw);
         tt_lmb{tabidx}.r = sum(rw)/den_r;
         % fusion of existence probabilities
     end
end