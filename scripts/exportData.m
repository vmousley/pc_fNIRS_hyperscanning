%% Export WTC Analysis Data
%{

Written by Victoria St. Clair, Centre for Brain and Cognitive Development
Email: v.stclair@bbk.ac.uk

This script exports WTC pipeline output for analysis in R.

Must specify line 22 depending on what kinds of data you want to export.
    - 'FINAL*_sscRegressed_hbo*' for hbo data with SSC regression
    - 'FINAL*_sscRegressed_hbr*' for hbr data with SSC regression
    - '*_raw_hbo*' for hbo data without SSC regression
    - '*_raw_hbr*' for hbr data without SSC regression
    - '*PDA*hbo*' for hbo pseudodyad data
    - '*PDA*hbr*' for hbo pseudodyad data

Change line 27 so file name matches type of export.

%}

clear all
clear
clc

config; 

exportfilePattern = '*FINAL_sscRegressed_hbo*'; % <<<< CHANGE depending on desired data

fullExportPath = fullfile(analysisPath, 'Hbo.mat');
% ^^^^^ CHANGE name depending on desired file name

fileList = dir(fullfile(analysisPath, exportfilePattern));
dataAll = [];

for i = 1:length(fileList)

    fileName = fileList(i).name;
    filePath = fullfile(exportPath, fileName);
    
    data = load(filePath);
    fields = fieldnames(data);

    dataAll.Dyads{i} = strrep(data.(fields{1}).Pnum,'-','');
    dataAll.Type{i} = data.(fields{1}).type;
    
    dataAll.Channel_ConditionMeans{i} = data.(fields{1}).Channel_ConditionMeans;
    dataAll.ROI_ConditionMeans{i} = data.(fields{1}).ROI_ConditionMeans;

end

save(fullExportPath, '-struct', 'dataAll');
