%reaction times 
function [go_licks_all,ngo_licks_all,dprimes,rt_go_all,rt_ngo_all,licks_raster,c_bias] = GNG_lick_rt_dprime (index_stim, lick_times, stim_times, stim_types, trial_responses, stim_ids, tone_dur , response_window,stim)


                    n_trials = 1:length(stim_times) ;
                    for t = 1:length(n_trials) ;
                        licks_after_stim = lick_times (lick_times > stim_times(t)) ;
                        licks_per_trial = licks_after_stim(licks_after_stim < stim_times(t)+tone_dur+response_window) ;
                        licks_from_stim_onset = licks_per_trial - stim_times(t) ;
                        licks_raster(t,1:length(licks_from_stim_onset)) = licks_from_stim_onset ;
                    end
                    licks_raster = licks_raster(:,1)';
                    
                    idx_go = find (stim_types(index_stim) ==  1);
                    idx_ngo = find (stim_types(index_stim) == -1);
                    
                    rt_go = mean (licks_raster(idx_go));
                    rt_ngo = mean (licks_raster(idx_ngo));

                    stim_types_level = stim_types(index_stim);
                    trial_responses_level = trial_responses(index_stim);
                    
                    binsize = length(stim_types_level)-1;
                    [~, ~, ~, ~, go_licks, ngo_licks] = GNG_trial_outcomes (binsize, stim_types_level, trial_responses_level);
                    [~, ~, dprime, c_bias] = GNG_dprime(binsize, go_licks, ngo_licks);
              
                    go_licks_all = mean(go_licks);
                    ngo_licks_all = mean(ngo_licks);
                    dprimes = mean(dprime);
                    rt_go_all = rt_go;
                    rt_ngo_all = rt_ngo;
end 