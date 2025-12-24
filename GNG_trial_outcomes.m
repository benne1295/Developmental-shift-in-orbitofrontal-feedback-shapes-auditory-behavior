function [hit, miss, fa, cr, go_licks, ngo_licks] = GNG_trial_outcomes (binsize, stim_types, trial_responses);

bining=zeros(1,length(trial_responses));
binind=1:binsize:length(bining);
for b=1:length(binind)-1;
    bining(binind(b):end)=bining(binind(b):end)+1;
end
bining=bining';


U = unique(bining,'stable');

for h= 1:length (U)
    hit(h)= sum (bining(1,:) == h & stim_types == 1 & trial_responses == 1);
    miss(h)=sum(bining(1,:)== h & stim_types == 1 & trial_responses == 0);
    fa(h)=sum(bining(1,:)== h & stim_types == -1 & trial_responses == 1);
    cr(h)=sum(bining(1,:)== h & stim_types == -1 & trial_responses == 0);
    go_licks(h) =  hit(h) / (hit(h)+ miss(h));
    ngo_licks(h) = fa(h) / (fa(h)+ cr(h));
    
    
if isnan(go_licks(h)) | isinf(go_licks(h))
    go_licks = 0 ;
end

if isnan(ngo_licks(h)) | isinf(ngo_licks(h))
    ngo_licks = 0 ;
end
    
    
end



end




