function [rasterX, rasterY] = GNG_raster (binArray, binSize, binCenters, rasterScale);

%rasterX = x vector of the rasterplot
%rasterY = y vector of the rasterplot 
[tr,b] = find(binArray); % find all spikes of the binarray
[rasterX,yy] = rasterize(binCenters(b)); 
rasterY = yy+reshape(repmat(tr',3,1),1,length(tr)*3); %reshape according to timebin of spike occurence

% scale the raster ticks
rasterY(2:3:end) = rasterY(2:3:end)+rasterScale; %scale raster Y according to smoothing 
end 


