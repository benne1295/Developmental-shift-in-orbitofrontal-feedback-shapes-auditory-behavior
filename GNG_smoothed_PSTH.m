function [fr_trial_smoothed, psth_Smoothed,binArray_smoothed,spikeCounts] = GNG_smoothed_PSTH (startRange, stopRange, binCenters, binArray, binSize, smoothSize);
    % psth_Smoothed = get the psth in hz
    % spikeCounts = get the spike sum per smoothed timebin 

    inclRange = binCenters > startRange & binCenters<= stopRange;
    spikeCounts = sum(binArray(:,inclRange),2)./(stopRange-startRange); % sum the spikecounts

    gaussian_window = gausswin(round(smoothSize*6),3); %apply a gaussian window 
    summed_window = gaussian_window./sum(gaussian_window); %sum all winows 
    binArray_smoothed = conv2(summed_window,1,binArray', 'same')'./binSize;  %conv the summed window 
    psth_Smoothed = mean(binArray_smoothed);
    bin_divided =  binArray_smoothed./binSize ;
    fr_trial_smoothed = mean (bin_divided') ;
    

end 