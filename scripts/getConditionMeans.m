%% Main WTC Analysis Script
%{

Written by Victoria Mousley, Centre for Brain and Cognitive Development
Email: v.mousley@bbk.ac.uk

This script averages WTC output across experimental conditions. 

Input: 
    data = analysisHbo or analysisHbr depending on which you want

%}

function getConditionMeans(data)

config; 

if isempty(data)
    disp('No data to input. Script will return.');
    return;
end

emptyChans = zeros(length(data.Channels_Rsq), width(data.Channels_Rsq));
meanCollabChans = repmat({emptyChans}, 1, length(data.End_Collaboration));
meanCollabScreenChans = repmat({emptyChans}, 1, length(data.End_CollaborationScreen));
meanIndividualChans = repmat({emptyChans}, 1, length(data.End_Individual));
meanRestChans = repmat({emptyChans}, 1, length(data.End_Rest));

emptyROIs = zeros(length(data.ROIs_Rsq), width(data.ROIs_Rsq));
meanCollabROIs = repmat({emptyROIs}, 1, length(data.End_Collaboration));
meanCollabScreenROIs = repmat({emptyROIs}, 1, length(data.End_CollaborationScreen));
meanIndividualROIs = repmat({emptyROIs}, 1, length(data.End_Individual));
meanRestROIs = repmat({emptyROIs}, 1, length(data.End_Rest));

%% For channels

for k = 1:numel(data.End_Collaboration)

    startIndex = data.Start_Collaboration(k);
    endIndex = startIndex + ((samp*cond_len)-1);

    for i=1:length(data.Channels_Rsq)
        for j = 1:width(data.Channels_Rsq)
            if all(isnan(data.Channels_Rsq{i,j}))
                meanCollabChans{1,k}(i,j) = NaN;
            else
                meanCollabChans{1,k}(i,j) = nanmean(data.Channels_Rsq{i,j}(:,startIndex:endIndex),2);
            end
        end
    end
end

for k = 1:numel(data.End_CollaborationScreen)

    startIndex = data.Start_CollaborationScreen(k);
    endIndex = startIndex + ((samp*cond_len)-1);

    for i=1:length(data.Channels_Rsq)
        for j = 1:width(data.Channels_Rsq)
            if all(isnan(data.Channels_Rsq{i,j}))
                meanCollabScreenChans{1,k}(i,j) = NaN;
            else
                meanCollabScreenChans{1,k}(i,j) = nanmean(data.Channels_Rsq{i,j}(:,startIndex:endIndex),2);
            end
        end
    end
end

for k = 1:numel(data.End_Individual)

    startIndex = data.Start_Individual(k);
    endIndex = startIndex + ((samp*cond_len)-1);

    for i=1:length(data.Channels_Rsq)
        for j = 1:width(data.Channels_Rsq)
            if all(isnan(data.Channels_Rsq{i,j}))
                meanIndividualChans{1,k}(i,j) = NaN;
            else
                meanIndividualChans{1,k}(i,j) = nanmean(data.Channels_Rsq{i,j}(:,startIndex:endIndex),2);
            end
        end
    end
end

for k = 1:numel(data.End_Rest)

    startIndex = data.Start_Rest(k);
    endIndex = startIndex + ((samp*rest_len)-1);

    for i=1:length(data.Channels_Rsq)
        for j = 1:width(data.Channels_Rsq)
            if all(isnan(data.Channels_Rsq{i,j}))
                meanRestChans{1,k}(i,j) = NaN;
            else
                meanRestChans{1,k}(i,j) = nanmean(data.Channels_Rsq{i,j}(:,startIndex:endIndex),2);
            end
        end
    end
end

% For region comparisons
for k = 1:numel(data.End_Collaboration)

    startIndex = data.Start_Collaboration(k);
    endIndex = startIndex + ((samp*cond_len)-1);

    for i=1:length(data.ROIs_Rsq)
        for j = 1:width(data.ROIs_Rsq)
            if all(isnan(data.ROIs_Rsq{i,j}))
                meanCollabROIs{1,k}(i,j) = NaN;
            elseif numel(data.ROIs_Rsq{i,j}) < endIndex
                meanCollabROIs{1,k}(i,j) = NaN;
            else
                meanCollabROIs{1,k}(i,j) = nanmean(data.ROIs_Rsq{i,j}(:,startIndex:endIndex),2);
            end
        end
    end
end

for k = 1:numel(data.End_CollaborationScreen)

    startIndex = data.Start_CollaborationScreen(k);
    endIndex = startIndex + ((samp*cond_len)-1);

    for i=1:length(data.ROIs_Rsq)
        for j = 1:width(data.ROIs_Rsq)
            if all(isnan(data.ROIs_Rsq{i,j}))
                meanCollabScreenROIs{1,k}(i,j) = NaN;
            elseif numel(data.ROIs_Rsq{i,j}) < endIndex
                meanCollabScreenROIs{1,k}(i,j) = NaN;
            else
                meanCollabScreenROIs{1,k}(i,j) = nanmean(data.ROIs_Rsq{i,j}(:,startIndex:endIndex),2);
            end
        end
    end
end

for k = 1:numel(data.End_Individual)

    startIndex = data.Start_Individual(k);
    endIndex = startIndex + ((samp*cond_len)-1);

    for i=1:length(data.ROIs_Rsq)
        for j = 1:width(data.ROIs_Rsq)
            if all(isnan(data.ROIs_Rsq{i,j}))
                meanIndividualROIs{1,k}(i,j) = NaN;
            elseif numel(data.ROIs_Rsq{i,j}) < endIndex
                meanIndividualROIs{1,k}(i,j) = NaN;
            else
                meanIndividualROIs{1,k}(i,j) = nanmean(data.ROIs_Rsq{i,j}(:,startIndex:endIndex),2);
            end
        end
    end
end

for k = 1:numel(data.End_Rest)

    startIndex = data.Start_Rest(k);
    endIndex = startIndex + ((samp*rest_len)-1);

    for i=1:length(data.ROIs_Rsq)
        for j = 1:width(data.ROIs_Rsq)
            if all(isnan(data.ROIs_Rsq{i,j}))
                meanRestROIs{1,k}(i,j) = NaN;
            elseif numel(data.ROIs_Rsq{i,j}) < endIndex
                meanRestROIs{1,k}(i,j) = NaN;
            else
                meanRestROIs{1,k}(i,j) = nanmean(data.ROIs_Rsq{i,j}(:,startIndex:endIndex));
            end
        end
    end
end

data.Channel_ConditionMeans = struct();
data.Channel_ConditionMeans.Collaboration = meanCollabChans;
data.Channel_ConditionMeans.CollaborationScreen = meanCollabScreenChans;
data.Channel_ConditionMeans.Individual = meanIndividualChans;
data.Channel_ConditionMeans.Rest = meanRestChans;

data.ROI_ConditionMeans = struct();
data.ROI_ConditionMeans.Collaboration = meanCollabROIs;
data.ROI_ConditionMeans.CollaborationScreen = meanCollabScreenROIs;
data.ROI_ConditionMeans.Individual = meanIndividualROIs;
data.ROI_ConditionMeans.Rest = meanRestROIs;

if strcmp(data.type, 'hbo')
    data.type = 'hbo';
    assignin('caller','finalHbo', data);
else
    strcmp(data.type, 'hbr');
    data.type = 'hbr';
    assignin('caller','finalHbr', data);
end

end