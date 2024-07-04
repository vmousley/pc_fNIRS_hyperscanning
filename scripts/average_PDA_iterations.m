% average WTC values obtained for set number iterations for PSEUDO dyads

function average_PDA_iterations(data, PDAmaxIterNum, temp, type)
        
    % initialize data
    PDAfinal = data;
    
    for trialNum=1:size(data.Channel_ConditionMeans.Collaboration,2) % loop over condition trials
            for iterNum = 1:PDAmaxIterNum % loop over iterations
                if iterNum==1
                    sum_CollabChans{trialNum}=temp.CollabChans{1,iterNum}{1,trialNum};
                    sum_CollabScreenChans{trialNum}=temp.CollabScreenChans{1,iterNum}{1,trialNum};
                    sum_IndividualChans{trialNum}=temp.IndividualChans{1,iterNum}{1,trialNum};
                    sum_CollabROIs{trialNum}=temp.CollabROIs{1,iterNum}{1,trialNum};
                    sum_CollabScreenROIs{trialNum}=temp.CollabScreenROIs{1,iterNum}{1,trialNum};
                    sum_IndividualROIs{trialNum}=temp.IndividualROIs{1,iterNum}{1,trialNum};
                else
                    sum_CollabChans{trialNum}=sum(cat(3,sum_CollabChans{trialNum},temp.CollabChans{1,iterNum}{1,trialNum}),3,'omitnan');
                    sum_CollabScreenChans{trialNum}=sum(cat(3,sum_CollabScreenChans{trialNum},temp.CollabScreenChans{1,iterNum}{1,trialNum}),3,'omitnan');
                    sum_IndividualChans{trialNum}=sum(cat(3,sum_IndividualChans{trialNum},temp.IndividualChans{1,iterNum}{1,trialNum}),3,'omitnan');
                    sum_CollabROIs{trialNum}=sum(cat(3,sum_CollabROIs{trialNum},temp.CollabROIs{1,iterNum}{1,trialNum}),3,'omitnan');
                    sum_CollabScreenROIs{trialNum}=sum(cat(3,sum_CollabScreenROIs{trialNum},temp.CollabScreenROIs{1,iterNum}{1,trialNum}),3,'omitnan');
                    sum_IndividualROIs{trialNum}=sum(cat(3,sum_IndividualROIs{trialNum},temp.IndividualROIs{1,iterNum}{1,trialNum}),3,'omitnan');
                end
            end
          RsqFinal_CollabChans{trialNum} = sum_CollabChans{trialNum}./iterNum;
          RsqFinal_CollabScreenChans{trialNum} = sum_CollabScreenChans{trialNum}./iterNum;
          RsqFinal_IndividualChans{trialNum} = sum_IndividualChans{trialNum}./iterNum;
          RsqFinal_CollabROIs{trialNum} = sum_CollabROIs{trialNum}./iterNum;
          RsqFinal_CollabScreenROIs{trialNum} = sum_CollabScreenROIs{trialNum}./iterNum;
          RsqFinal_IndividualROIs{trialNum} = sum_IndividualROIs{trialNum}./iterNum;
    end
     
    for restNum=1:size(data.Channel_ConditionMeans.Rest,2)
        for iterNum = 1:PDAmaxIterNum
            if iterNum==1
                sum_RestChans{restNum}=temp.RestChans{1,iterNum}{1,restNum};
                sum_RestROIs{restNum}=temp.RestROIs{1,iterNum}{1,restNum};
            else
                sum_RestChans{restNum}=sum(cat(3,sum_RestChans{restNum},temp.RestChans{1,iterNum}{1,restNum}),3,'omitnan');
                sum_RestROIs{restNum}=sum(cat(3,sum_RestROIs{restNum},temp.RestROIs{1,iterNum}{1,restNum}),3,'omitnan');
            end
        end
      RsqFinal_RestChans{restNum} = sum_RestChans{restNum}./iterNum;
      RsqFinal_RestROIs{restNum} = sum_RestROIs{restNum}./iterNum;
    end

    result.info = ['This struct includes pseudo-dyad wavelet transform ' ...
    'coherence (WTC) calculations for every possible channel and ROI ' ...
    'combination. WTC was calculated between p1 and p2 after the phase ' ...
    'scrambling of p2 signals. Data were averaged over 100 iterations.' ...
    'For child-mum study, p1 = child and p2 = mum. RsqFinal houses WTC ' ...
    'for each p1 (rows) and p2 (columns) pairing. Excluded channels ' ...
    'are identified via an Excel sheet of manual identification and ' ...
    'include SSCs to be removed. For these channels, RsqFinal are' ...
    'NaNs. RsqFinal values are for haemoglobin type in .type and' ...
    'are averaged across the period of interest (0.02 to 0.10Hz). Email' ...
    'v.mousley@bbk.ac.uk with queries.'];

    PDAfinal.Channel_ConditionMeans.Collaboration = RsqFinal_CollabChans;
    PDAfinal.Channel_ConditionMeans.CollaborationScreen = RsqFinal_CollabScreenChans;
    PDAfinal.Channel_ConditionMeans.Individual = RsqFinal_IndividualChans;
    PDAfinal.Channel_ConditionMeans.Rest = RsqFinal_RestChans;
    PDAfinal.ROI_ConditionMeans.Collaboration = RsqFinal_CollabROIs;
    PDAfinal.ROI_ConditionMeans.CollaborationScreen = RsqFinal_CollabScreenROIs;
    PDAfinal.ROI_ConditionMeans.Individual = RsqFinal_IndividualROIs;
    PDAfinal.ROI_ConditionMeans.Rest = RsqFinal_RestROIs;

    if strcmp(type, 'hbo')
        PDAfinal.type = 'hbo';
        assignin('base','PDAfinalHbo', PDAfinal);
    else
        strcmp(type, 'hbr');
        PDAfinal.type = 'hbr';
        assignin('base','PDAfinalHbr', PDAfinal);
    end
end