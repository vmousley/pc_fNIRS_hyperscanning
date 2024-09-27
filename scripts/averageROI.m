% Average Within-ROI
%{

Written by Victoria St. Clair, Centre for Brain and Cognitive Development
Email: v.stclair@bbk.ac.uk

This function averages available long-separation channels within ROIs.

Inputs: 
    data = fNIRS data
    region = region var as specified in config.m
    regionName = string for region name

%}

function averageROI (data, region, regionName)

config;

concatDoublesP1 = [];
concatDoublesP2 = [];

% p1 averaging by region
for i = 1:numel(data.channels)
    if ismember(i, region)
        if ismember(i, data.p1ExcludedChannels)
            doubleArrayP1 = NaN(size(data.rawP1(:,i)));
            concatDoublesP1 = [concatDoublesP1, doubleArrayP1];
        else
            doubleArrayP1 = data.rawP1(:,i);
            concatDoublesP1 = [concatDoublesP1, doubleArrayP1];
        end
    end
end

if sum(any(~isnan(concatDoublesP1)))<2
    concatDoublesP1 = NaN(length(data.rawP1),1);
    p1MeanROI = nanmean(concatDoublesP1,2);
else
    p1MeanROI = nanmean(concatDoublesP1,2);
end

% p2 averaging by region
for i = 1:numel(data.channels)
    if ismember(i, region)
        if ismember(i, data.p2ExcludedChannels)
            doubleArrayP2 = NaN(size(data.rawP2(:,i)));
            concatDoublesP2 = [concatDoublesP2, doubleArrayP2];
        else
            doubleArrayP2 = data.rawP2(:,i);
            concatDoublesP2 = [concatDoublesP2, doubleArrayP2];
        end
    end
end

if sum(any(~isnan(concatDoublesP2)))<2
    concatDoublesP2 = NaN(length(data.rawP2),1);
    p2MeanROI = nanmean(concatDoublesP2,2);
else
    p2MeanROI = nanmean(concatDoublesP2,2);
end

newNameP1 = genvarname([regionName '_p1']);
newNameP2 = genvarname([regionName '_p2']);

assignin('caller',newNameP1, p1MeanROI)
assignin('caller',newNameP2, p2MeanROI)

end
