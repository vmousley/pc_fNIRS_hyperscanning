%% Main WTC Script
%{

Written by Victoria Mousley, Centre for Brain and Cognitive Development
and Letizia Contini, Politecnico di Milano
Email: v.mousley@bbk.ac.uk

This script runs the main fNIRS WTC pipeline supporting Mousley et al. 
(under review). Analytical pipeline optimisation in developmental fNIRS 
hyperscanning data: Neural coherence between preschoolers collaborating 
with their mothers. 

Input = pre-processed fNIRS data from children and mothers
Output = true and phase-scrambled pseudodyad WTC data as labelled .structs

Must specify:
    - Line 31 for number of desired pseudo dyad iterations. If not running
    phase-scrambled pseudodyad analysis, set to 0.
    - Line 32 for SSC regression if desired ('true' if yes, 'false' if
    not). NB: Setting ssc reg = true means physiological regression will
    be used where dyads have at least one valid SSC. Otherwise, unregressed
    data will be fed into WTC and results will be saved as 'raw'.

%}

close all
clear 
clc

%% SET DESIRED PARAMETERS %%
totalIterations=0; % <<<< SET to number of desired phase-scrambling iterations
SSCR_ON = true; % <<<< SET to true if SSC regression desired, otherwise false
 
%%
config;
inputPath = basePath;
outputPath = analysisPath;
dyadFolders = dir(fullfile(inputPath, 'Dyad*'));
totalDyadsNumber=size(dyadFolders,1);

if SSCR_ON == true
    sscReg = 'sscRegressed';
else
    sscReg = 'raw';
end

%%
for i = 1:numel(dyadFolders) 
    dyadFolder = dyadFolders(i).name;
    dyadFolderPath = fullfile(inputPath, dyadFolder);
    cFolder = fullfile(dyadFolderPath, sprintf('C%s', dyadFolder(5:end))); % child data
    mFolder = fullfile(dyadFolderPath, sprintf('M%s', dyadFolder(5:end))); % mother data
    cFiles = dir(fullfile(cFolder, 'C*_ppr.mat'));
    mFiles = dir(fullfile(mFolder, 'M*_ppr.mat'));

    % check for data belonging to the same dyad
    if numel(cFiles) == 1 && numel(mFiles) == 1  

        % load data
        data1 = load(fullfile(cFolder, cFiles.name));
        data2 = load(fullfile(mFolder, mFiles.name));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% calculate TRUE dyads WTC %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % calculate channel-level WTC for Hbo
        channelWTC(data1, data2, 'hbo', SSCR_ON)

        % if ssc reg is true but no valid sscs, reset to false
        if SSCR_ON == 1
            if isempty(resultHbo.p1SSlistGood) && isempty(resultHbo.p2SSlistGood)
                SSCR_ON = 0;
                sscReg = 'raw';
                fprintf('SSC reg set but no valid data. Resetting SSC arg to false. \n')
            else
                fprintf('Proceeding into pipeline with SSC regression. \n')
            end
        elseif SSCR_ON == 0
            fprintf('SSC regression false. Proceeding without regression. \n')
        end

        % calculate ROI-wise WTC
        roiWTC(resultHbo, SSCR_ON);

        % average channel- and ROI-wise data across conditions
        getConditionMeans(analysisHbo);

        % save data 
        dyadName = analysisHbo.Pnum;
        fileNameHbo = sprintf('FINAL_%s_%s_hbo.mat', dyadName, sscReg);
        save(fullfile(outputPath, fileNameHbo), 'finalHbo');
        
        % repeat for Hbr
        channelWTC(data1, data2, 'hbr', SSCR_ON);
        roiWTC(resultHbr, SSCR_ON);
        getConditionMeans(analysisHbr);
        dyadName = analysisHbr.Pnum;
        fileNameHbr = sprintf('FINAL_%s_%s_hbr.mat', dyadName, sscReg);
        save(fullfile(outputPath, fileNameHbr), 'finalHbr');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%% calculate PSEUDO dyads WTC %%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if totalIterations > 0
            startTime = tic; % to monitor computation time
            tempHbo = struct();
            tempHbr = struct();

            for iterNum = 1:totalIterations % loop over phase-scrambling iterations
            
                iterStartTime = tic; % to monitor computation time (per iteration)
            
                % phase-scrambling of mother data
                data2_scrambled = phase_scrambling_dc(data2);
            
                % calculate channel-level WTC for Hbo
                clear resultHbo; % clear results from the previous iteration
                channelWTC(data1, data2_scrambled, 'hbo', SSCR_ON)

                if exist('resultHbo', 'var')
                    roiWTC(resultHbo, SSCR_ON);
                    getConditionMeans(analysisHbo);
                    dyadName = analysisHbo.Pnum;    
                else
                    fprintf('%s does not have any valid channels for one or both participants. Script moving on. \n', dyadFolder);
                end 

                % save temporal data - needed to extract average values
                tempHbo.CollabChans{iterNum} = finalHbo.Channel_ConditionMeans.Collaboration;
                tempHbo.CollabScreenChans{iterNum} = finalHbo.Channel_ConditionMeans.CollaborationScreen;
                tempHbo.IndividualChans{iterNum} = finalHbo.Channel_ConditionMeans.Individual;
                tempHbo.RestChans{iterNum} = finalHbo.Channel_ConditionMeans.Rest;
                tempHbo.CollabROIs{iterNum} = finalHbo.ROI_ConditionMeans.Collaboration;
                tempHbo.CollabScreenROIs{iterNum} = finalHbo.ROI_ConditionMeans.CollaborationScreen;
                tempHbo.IndividualROIs{iterNum} = finalHbo.ROI_ConditionMeans.Individual;
                tempHbo.RestROIs{iterNum} = finalHbo.ROI_ConditionMeans.Rest;
    
                % calculate Channel-level WTC for Hbr
                clear resultHbr; % clear results from the previous iteration
                channelWTC(data1,data2_scrambled, 'hbr',SSCR_ON)
    
                if exist('resultHbr', 'var')
                    % calculate ROI-level WTC for Hbr
                    roiWTC(resultHbr,SSCR_ON);
                    % extract mean WTC for each condition
                    getConditionMeans(analysisHbr);
                    dyadName = analysisHbr.Pnum;    
                else
                    fprintf('%s does not have any valid channels for one or both participants. Script moving on. \n', dyadFolder);
                end
                % save temporal data - needed to extract average values
                tempHbr.CollabChans{iterNum} = finalHbr.Channel_ConditionMeans.Collaboration;
                tempHbr.CollabScreenChans{iterNum} = finalHbr.Channel_ConditionMeans.CollaborationScreen;
                tempHbr.IndividualChans{iterNum} = finalHbr.Channel_ConditionMeans.Individual;
                tempHbr.RestChans{iterNum} = finalHbr.Channel_ConditionMeans.Rest;
                tempHbr.CollabROIs{iterNum} = finalHbr.ROI_ConditionMeans.Collaboration;
                tempHbr.CollabScreenROIs{iterNum} = finalHbr.ROI_ConditionMeans.CollaborationScreen;
                tempHbr.IndividualROIs{iterNum} = finalHbr.ROI_ConditionMeans.Individual;
                tempHbr.RestROIs{iterNum} = finalHbr.ROI_ConditionMeans.Rest;
    
                % estimate remaining time
                elapsedTime = toc(iterStartTime);
                estimatedTimeRemaining = ((totalIterations-iterNum) * elapsedTime) / 60;
    
                % show progress and estimated remaining time for PSEUDO dyad analysis
                fprintf('Dyad %i/%i. Iteration %i/%i completed. Estimated %.2f minutes remaining. \n\n', ...
                     i, totalDyadsNumber, iterNum, totalIterations, estimatedTimeRemaining);
            end
            TotCompTime = toc(startTime); % save total computation time
    
            % average PSEUDO dyads WTC over x iterations 
            average_PDA_iterations(finalHbo, totalIterations, tempHbo, 'hbo');
            fileNameHbo = sprintf('PDA_%s_%s_hbo.mat', PDAfinalHbo.Pnum, sscReg);
                save(fullfile(outputPath, fileNameHbo), 'PDAfinalHbo');
    
            average_PDA_iterations(finalHbr, totalIterations, tempHbr, 'hbr');
            fileNameHbr = sprintf('PDA_%s_%s_hbr.mat', PDAfinalHbr.Pnum, sscReg);
                save(fullfile(outputPath, fileNameHbr), 'PDAfinalHbr');
    
            % clear temporary variables
            clear tempHbo
            clear tempHbr
    
            % estimate remaining time
            estimatedDaysRemaining = ((totalDyadsNumber-i) * TotCompTime) / (60*60*24);
    
            % show progress and estimated remaining time
            fprintf('\n###### %s completed! ###### \n', dyadFolder)
            fprintf('Dyads completed: %i/%i. Estimated %.2f days remaining. \n\n', i, totalDyadsNumber, estimatedDaysRemaining)
        else
            fprintf('Skipping pseudodyad analysis. \n')
        end

    else
        fprintf('Missing or multiple files in C and M folders. Skipping %s folder. \n', dyadFolder);
    end

end
