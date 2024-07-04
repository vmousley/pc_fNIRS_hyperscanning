% The SD variable is the SD structure that you should already have
% Put HbO and HbR in the right format for the dc variable;
% HbO is [samples x channels]

dc(:,:,1)=HbO;
dc(:,:,2)=HbR;

% Write here the numbers of the good SSCs
SSlistGood=[3 16];

% filtered_dc is the concentration data after the regression of SSCs
[filtered_dc,Stats] = PhysiologyRegression_GLM_fnirs_course(dc,SD,SSlistGood,[]);

%%%% squeeze(filtered_dc(:,:,1) is the HbO to give as an input to WTC
%%%% squeeze(filtered_dc(:,:,2) is the HbR to give as an input to WTC