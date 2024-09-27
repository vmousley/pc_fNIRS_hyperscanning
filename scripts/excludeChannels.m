%% Exclude Noisy Channels
%{

Written by Victoria St. Clair, Centre for Brain and Cognitive Development
Email: v.stclair@bbk.ac.uk

This script identifies noisy channels specified in Excel sheet
to exclude before WTC calculation for expediency. 

Inputs: 
    data = NIRS data
    varName = variable name for data

%}

function excludeChannels(data, varName)

config;

filename ='ChannelExclusions.xlsx';

if contains(data.Pnum,'C')
    sheet = 'data1_exclusions';
else 
    contains(data.Pnum, 'M');
    sheet = 'data2_exclusions';
end

numericPart = regexp(data.Pnum, '\d+', 'match');
numericValue = str2double(numericPart);

channels = readtable(filename, 'Sheet', sheet,'ReadVariableNames',true);
idx = find(strcmp(channels.Properties.VariableNames,data.Pnum));
excludedChannels = find(channels{:,(idx)}==0);

SSlistGood = setdiff(sscs, excludedChannels);

data.allExcludedChannels = sortrows(unique(vertcat(excludedChannels, sscs)));
data.SSlistGood = SSlistGood;

assignin('caller', varName, data);

end
