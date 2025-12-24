function [fr_trial, psth, binArray, binCenters] = GNG_binning_PSTH (startRange,stopRange, binSize, eventTimes, spikeTimes);

%psth = peri stimulus time histogram of the chosen unit
%binArray = binary array of all spikes per trial 
%binCenters = binary output 

binBorders = startRange:binSize:stopRange ; %define the time points 
numBins = length(binBorders)-1 ; 
binArray = zeros(length(eventTimes), numBins) ; 

for r = 1:length(eventTimes)
    [n,binCenters] = histdiff(spikeTimes, eventTimes(r), binBorders); %get the spike per binsize
    binArray(r,:) = n;
end

  psth = mean(binArray./binSize) ; % normalize to Hz           
  fr_trial =  binArray./binSize ;
  fr_trial = mean(fr_trial') ;
  
end
