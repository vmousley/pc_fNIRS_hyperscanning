% ROI-Wise WTC
%{

Written by Victoria Mousley, Centre for Brain and Cognitive Development
Email: v.mousley@bbk.ac.uk

This script averages fNIRS signals within regions across available 
long-separation channels.

Inputs: 
    data = 'resultHbo' or 'resultHbr', depending on which you want
    ssc_reg = 'true' if running analysis with regression, 'false' if without
%}

function roiWTC(data, ssc_reg)

config;

if isempty(data)
    disp('No data to input. Script will return.');
    return;
end

averageROI(data, leftPFC, 'leftPFCMean')
averageROI(data, rightPFC, 'rightPFCMean')
averageROI(data, leftTPJ, 'leftTPJMean')
averageROI(data, rightTPJ, 'rightTPJMean')

withinROI.labels = cell(1,4);
withinROI.p1 = cell(1,4);
withinROI.p2 = cell(1,4);

for i = 1:4
    if i == 1
        withinROI.p1{i} = leftPFCMean_p1;
        withinROI.p2{i} = leftPFCMean_p2;
        withinROI.labels{i} = 'leftPFC';
    elseif i == 2
        withinROI.p1{i} = rightPFCMean_p1;
        withinROI.p2{i} = rightPFCMean_p2;
        withinROI.labels{i} = 'rightPFC';
    elseif i == 3
        withinROI.p1{i} = leftTPJMean_p1;
        withinROI.p2{i} = leftTPJMean_p2;
        withinROI.labels{i} = 'leftTPJ';
    elseif i == 4
        withinROI.p1{i} = rightTPJMean_p1;
        withinROI.p2{i} = rightTPJMean_p2;
        withinROI.labels{i} = 'rightTPJ';
    else 
        print('ERROR!')
    end
end

p1NumROIs = width(withinROI.p1);
p2NumROIs = width(withinROI.p2);

RsqFinal = cell(p1NumROIs, p2NumROIs);
pnoi = zeros(2,1);

p1ROIs = false(size(withinROI.p1));
for i = 1:numel(withinROI.p1)
    if ~any(isnan(withinROI.p1{i}))
        p1ROIs(i) = true;
    end
end

p2ROIs = false(size(withinROI.p2));
for k = 1:numel(withinROI.p2)
    if ~any(isnan(withinROI.p2{k}))
        p2ROIs(k) = true;
    end
end

for col1 = 1:p1NumROIs
    for col2 = 1:p2NumROIs

        if any(isnan(withinROI.p1{:,col1}))
            RsqFinal{col1,col2} = NaN(size(data.t));
        elseif any(isnan(withinROI.p2{:,col2}))
            RsqFinal{col1,col2} = NaN(size(data.t));
        else
            data1_wtc = [data.t', withinROI.p1{col1}];
            data2_wtc = [data.t', withinROI.p2{col2}];

            [Rsq, period, ~, ~, ~] = wtc(data1_wtc, data2_wtc, 'mcc', 0);

            pnoi(1) = find(period > poi(1), 1, 'first');
            pnoi(2) = find(period < poi(2), 1, 'last');

            RsqMeanCoherence = mean(Rsq(pnoi(1):pnoi(2),:));

            RsqFinal{col1,col2} = RsqMeanCoherence;
            
        end
    end
end

data.ROIs_Rsq = RsqFinal;
data.ROI_labels = withinROI.labels;
data.ROIs_info = ['ROIs analysis has p1 ROIs as rows and p2s as columns. ' ...
    'ROIs according to labels field.'];

if strcmp(data.type, 'hbo')
    data.type = 'hbo';
    assignin('caller','analysisHbo', data);
else
    strcmp(data.type, 'hbr');
    data.type = 'hbr';
    assignin('caller','analysisHbr', data);

end