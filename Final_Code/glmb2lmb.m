function tt_lmb= glmb2lmb(glmb)

%find unique labels (with different possibly different association histories)
lmat= zeros(2,length(glmb.tt),1); 
for tabidx= 1:length(glmb.tt)
    lmat(:,tabidx)= glmb.tt{tabidx}.l;
end
lmat= lmat'; %matrix of labels, each label is a row
%sort(lmat);

[cu,~,ic]= unique(lmat,'rows'); cu= cu'; %cu is a matrix, each label is a column

%initialize LMB struct
tt_lmb= cell(size(cu,2),1);

for tabidx=1:length(tt_lmb)
   tt_lmb{tabidx}.r= 0;
   tt_lmb{tabidx}.x= [];
   tt_lmb{tabidx}.w= []; %particle weights
   tt_lmb{tabidx}.l= cu(:,tabidx);
end

Num_particles = zeros(1,length(tt_lmb));
Particle_Index = ones(1,length(tt_lmb));

for hidx=1:length(glmb.w) %for each hypothesis  
   for t= 1:glmb.n(hidx)    %for each element in the hypothesis
      trkidx= glmb.I{hidx}(t);  %find the t'th track
      newidx= ic(trkidx);
      Num_particles(newidx)  = Num_particles(newidx) + length(glmb.tt{trkidx}.w); %num particles in each new lmb component
   end
end

for tabidx=1:length(tt_lmb)
   tt_lmb{tabidx}.x= zeros(Num_particles(tabidx),4);
   tt_lmb{tabidx}.w= zeros(Num_particles(tabidx),1);
end


%extract individual tracks
for hidx=1:length(glmb.w) %each hypothesis
   for t= 1:glmb.n(hidx) %each element in a hypothesis
      trkidx= glmb.I{hidx}(t); %element t of hypothesis 
      newidx= ic(trkidx);
            
      tt_lmb{newidx}.x(Particle_Index(newidx):(Particle_Index(newidx)+length(glmb.tt{trkidx}.w)-1),:) = glmb.tt{trkidx}.x;
      tt_lmb{newidx}.w(Particle_Index(newidx):(Particle_Index(newidx)+length(glmb.tt{trkidx}.w)-1)) ...
          = glmb.w(hidx)*glmb.tt{trkidx}.w;
      Particle_Index(newidx) = Particle_Index(newidx) + length(glmb.tt{trkidx}.w);
      %tt_lmb{newidx}.x= cat(2,tt_lmb{newidx}.x,glmb.tt{trkidx}.x);
      %tt_lmb{newidx}.w= cat(1,tt_lmb{newidx}.w,glmb.w(hidx)*glmb.tt{trkidx}.w);
   end
end

%extract existence probabilities and normalize track weights
for tabidx=1:length(tt_lmb)
   tt_lmb{tabidx}.r= sum(tt_lmb{tabidx}.w); %added 7/10/21. ensures no potential matlab rounding error to 1
   tt_lmb{tabidx}.w= tt_lmb{tabidx}.w/tt_lmb{tabidx}.r;
end

end