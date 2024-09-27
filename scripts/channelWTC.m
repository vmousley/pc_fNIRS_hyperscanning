%% Channel-Wise WTC
%{

Written by Victoria St. Clair, Centre for Brain and Cognitive Development
Email: v.stclair@bbk.ac.uk

This script takes estimates the wavelet transform coherence, calculated using
the analytic Morlet wavelet, between every possible valid channel-to-channel
pairing of two participants.

Inputs:

    data1 = fNIRS data from participant 1 (structured as in example_data files) 
    data2 = fNIRS data from participant 2 (structured as in example_data files)
    type = 'hbo' for oxygenated haemoglobin, 'hbr' for deoxygenated haemoglobin
    ssc_reg = 'true' if running analysis with regression, 'false' if without

%}

function channelWTC(data1, data2, type, ssc_reg)

config;

fprintf('Starting WTC function... \n');

if nargin < 4
    error('Not enough input arguments. \n');
end

validOptions = {'hbo', 'hbr'};

if ~ismember(type, validOptions)
    error('Invalid value for type. Valid options are: hbo or hbr. \n');
end

if isequal(data1.t, data2.t)
    fprintf('CHECK PASSED: Time vectors match. \n')
else
    fprintf('WARNING: timelines differ. \n')
end

if isequal(data1.Start_Collaboration, data2.Start_Collaboration) && ...
    isequal(data1.End_Collaboration, data2.End_Collaboration) && ...
    isequal(data1.Duration_Collaboration, data2.Duration_Collaboration)
    fprintf('CHECK PASSED: Collaboration condition event markers match. \n');
else 
    fprintf('WARNING: Collaboration condition event markers DO NOT match. \n');
end

if isequal(data1.Start_CollaborationScreen, data2.Start_CollaborationScreen) && ...
    isequal(data1.End_CollaborationScreen, data2.End_CollaborationScreen) && ...
    isequal(data1.Duration_CollaborationScreen, data2.Duration_CollaborationScreen)
    fprintf('CHECK PASSED: Collaboration Screen condition event markers match. \n');
else 
    fprintf('WARNING: Collaboration Screen condition event markers DO NOT match. \n');
end

if isequal(data1.Start_Individual, data2.Start_Individual) && ...
    isequal (data1.End_Individual, data2.End_Individual) && ...
    isequal(data1.Duration_Individual, data2.Duration_Individual)
    fprintf('CHECK PASSED: Individual condition event markers match. \n');
else 
    fprintf('WARNING: Individual condition event markers DO NOT match. \n');
end

if isequal(data1.Start_Rest, data2.Start_Rest) && ...
    isequal (data1.End_Rest, data2.End_Rest) && ...
    isequal(data1.Duration_Rest, data2.Duration_Rest)
    fprintf('CHECK PASSED: Rest condition event markers match. \n');
else 
    fprintf('WARNING: Rest condition event markers DO NOT match. \n')
end

excludeChannels(data1, 'data1');
excludeChannels(data2, 'data2');

if isempty(data1.SSlistGood) && isempty(data2.SSlistGood)
    ssc_reg = false;
    SSCR_ON = 0;
end

p1NumChans = data1.Ch;
p2NumChans = data2.Ch;

RsqFinal = cell(p1NumChans, p2NumChans);
pnoi = zeros(2,1);

%%
if ssc_reg == true
    dc_orig1 = data1.procResult.dc(:,1:2,:);
    dc1 = permute(dc_orig1, [1, 3, 2]); 

    dc_orig2 = data2.procResult.dc(:,1:2,:);
    dc2 = permute(dc_orig2, [1, 3, 2]);

    if isempty(data1.SSlistGood) && isempty(data2.SSlistGood)
        fprintf('SSC reg is true but neither participant has a valid channel. Continuing with raw data. \n');
        return;
    elseif isempty(data1.SSlistGood) && ~isempty(data2.SSlistGood)
        fprintf('SSC reg is true, but only data 2 has valid SSC. Running reg on data 2 only. \n');
        filtered_dc1 = dc1;
        [filtered_dc2] = PhysiologyRegression_GLM_fnirs_course(dc2, data2.SD, data2.SSlistGood,[]);
    elseif ~isempty(data1.SSlistGood) && isempty(data2.SSlistGood)
        fprintf('SSC reg is true, but only data 1 has valid SSC. Running reg on data 1 only. \n');
        [filtered_dc1] = PhysiologyRegression_GLM_fnirs_course(dc1, data1.SD, data1.SSlistGood,[]);
        filtered_dc2 = dc2;
    elseif ~isempty(data1.SSlistGood) && ~isempty(data2.SSlistGood)
        fprintf('SSC reg is true and both participants have at least one valid SSC. Running reg on both data 1 and data 2. \n');
        [filtered_dc1] = PhysiologyRegression_GLM_fnirs_course(dc1, data1.SD, data1.SSlistGood,[]);
        [filtered_dc2] = PhysiologyRegression_GLM_fnirs_course(dc2, data2.SD, data2.SSlistGood,[]);
    end

    if strcmp(type, 'hbo')
        p1AllChannels = squeeze(filtered_dc1(:,:,1));
        p2AllChannels = squeeze(filtered_dc2(:,:,1));

    elseif strcmp(type,'hbr')
        p1AllChannels = squeeze(filtered_dc1(:,:,2));
        p2AllChannels = squeeze(filtered_dc2(:,:,2));
    end

elseif ssc_reg == false
    disp('SSC regression is false, either because it was set 0 or because neither participant has a valid SSC data. Processing raw signal.')

    if strcmp(type, 'hbo')
        p1AllChannels = squeeze(data1.procResult.dc(:,1,:));
        p2AllChannels = squeeze(data2.procResult.dc(:,1,:));
    elseif strcmp(type, 'hbr')
        p1AllChannels = squeeze(data1.procResult.dc(:,2,:));
        p2AllChannels = squeeze(data2.procResult.dc(:,2,:));
    end

else
    error('ERROR: Invalid value for SSC_reg argument. Valid options are: true or false. \n');
end

%% Account for small discrepencies in time intervals

data1T=diff(data1.t(:,1));
data2T=diff(data2.t(:,1));

if any(abs(data1T-data1T(1))>1e-1*data1T(1))
    error('Time step in data 1 time must be constant. \n');
end

if any(abs(data2T-data2T(1))>1e-1*data2T(1))
    error('Time step in data 2 time must be constant. \n');
end

minData1 = min(data1.t);
maxData1 = max(data1.t);

minData2 = min(data2.t);
maxData2 = max(data2.t);

data1.timeRevised = (minData1:data1T:maxData1)';
data2.timeRevised = (minData2:data2T:maxData2)';

if length(data1.timeRevised) ~= length(data1.t)
    numRowsRemoved1 = abs(length(data1.timeRevised)-length(data1.t));
    numActualRows1 = length(p1AllChannels);
    newNumRows1 = numActualRows1 - numRowsRemoved1;
    p1AllChannelsRevised = p1AllChannels(1:newNumRows1,:);

    lastIndex = data1.End_Rest(end);
    data1.End_Rest(end) = lastIndex - numRowsRemoved1;

    fprintf('Revised timestamps of data1 do not match original. %.0f row(s) removed and rest index altered. \n', numRowsRemoved1);
else
    p1AllChannelsRevised = p1AllChannels(:,:);
    disp('Revised timestamps of data1 match original.')
end

if length(data2.timeRevised) ~= length(data2.t)
    numRowsRemoved2 = abs(length(data2.timeRevised)-length(data2.t));
    numActualRows2 = length(p2AllChannels);
    newNumRows2 = numActualRows2 - numRowsRemoved2;
    p2AllChannelsRevised = p2AllChannels(1:newNumRows2,:);

    lastIndex = data2.End_Rest(end);
    data2.End_Rest(end) = lastIndex - numRowsRemoved2;

    fprintf('Revised timestamps of data2 do not match original. %.0f row(s) removed and rest index altered. \n', numRowsRemoved2);
else
    p2AllChannelsRevised = p2AllChannels(:,:);
    disp('Revised timestamps of data2 match original.')
end

%% Calculate WTC

for col1 = 1:p1NumChans
    for col2 = 1:p2NumChans

        if nnz(ismember(data1.allExcludedChannels, col1))>0
            RsqFinal{col1,col2} = NaN(size(p1AllChannelsRevised(:,col1)'));
        elseif nnz(ismember(data2.allExcludedChannels, col2))>0
            RsqFinal{col1,col2} = NaN(size(p2AllChannelsRevised(:,col2)'));
        else
            data1_wtc = [data1.timeRevised, p1AllChannelsRevised(:, col1)];
            data2_wtc = [data2.timeRevised, p2AllChannelsRevised(:, col2)];

            [Rsq, period, ~, ~, ~] = wtc(data1_wtc, data2_wtc, 'mcc', 0);

            pnoi(1) = find(period > poi(1), 1, 'first'); %
            pnoi(2) = find(period < poi(2), 1, 'last');

            RsqMeanCoherence = mean(Rsq(pnoi(1):pnoi(2),:));

            RsqFinal{col1, col2} = RsqMeanCoherence;

        end
    end
end

result.info = ['This struct holds wavelet transform coherence (WTC) ' ...
    'calculations for every possible channel combination between p1 and p2. ' ...
    'For child-mum study, p1 = child and p2 = mum. RsqFinal = WTC ' ...
    'for each p1 (rows) and p2 (columns) pairing. Excluded channels ' ...
    'are identified via an Excel sheet of manual identification and ' ...
    'include SSCs to be removed. For these channels, RsqFinal are' ...
    'NaNs. RsqFinal values are for haemoglobin type in .type and' ...
    'are averaged across the period of interest specified in config. Email' ...
    'v.mousley@bbk.ac.uk with queries.'];

result.Pnum = [data1.Pnum, '-' data2.Pnum];

result.Channels_Rsq = RsqFinal;
result.Channels_Info = ['Channel Rsqs are p1 as rows and p2 as columns.'];
result.period = period;
result.pnoi = pnoi;
result.poi = poi;

result.channels = 1:1:size(RsqFinal,2);
result.p1ExcludedChannels = data1.allExcludedChannels;
result.p2ExcludedChannels = data2.allExcludedChannels;
result.t = data1.t';
result.rawP1 = p1AllChannels;
result.rawP2 = p2AllChannels;
result.ssc_reg = ssc_reg;
result.p1SSlistGood = data1.SSlistGood;
result.p2SSlistGood = data2.SSlistGood;

% this assumes that the condition timestamps are the same for p1 and p2, 
% which has been checked in the beginning 

result.CondNames = [data1.CondNames];
result.Duration_Collaboration = data1.Duration_Collaboration;
result.Duration_CollaborationScreen = data1.Duration_CollaborationScreen;
result.Duration_Individual = data1.Duration_Individual;
result.Duration_Rest = data1.Duration_Rest;

result.Start_Collaboration = data1.Start_Collaboration;
result.Start_CollaborationScreen = data1.Start_CollaborationScreen;
result.Start_Individual = data1.Start_Individual;
result.Start_Rest = data1.Start_Rest;

result.End_Collaboration = data1.End_Collaboration;
result.End_CollaborationScreen = data1.End_CollaborationScreen;
result.End_Individual = data1.End_Individual;
result.End_Rest = data1.End_Rest;

if strcmp(type, 'hbo')
    result.type = 'hbo';
    assignin('base','resultHbo', result);
else
    strcmp(type, 'hbr');
    result.type = 'hbr';
    assignin('base','resultHbr', result);
end
