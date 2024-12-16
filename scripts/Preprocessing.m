%% This script runs the fNIRS data pre-processing, and works on the .nirs file
%% It loads a .NIRS file, calculates the sampling frequency, and checks probe design
clear all
close all

%% Select .nirs file to preprocess

[filename_nirs, pathname] = uigetfile('*.nirs','Select the fNIRS data cut file');

% Extract participant number from file name
Pnum=filename_nirs(1:3);

% Create a full path to the file
filename_nirs=fullfile(pathname ,filename_nirs);

% Load .NIRS file
load(filename_nirs,'-mat')

% Calculate sampling frequency based on time stamps of the data
Data.fs=1/(t(2)-t(1));
1/(t(2)-t(1))

% Create list for number of channels per participant (18)
SD.MeasListAct=ones(18,1);

% Calculate the number of channels based on size of data matrix and the # of wavelengths
Data.nCh=size(d,2)/length(SD.Lambda); % Number of channels

%% Check probe design

% Calculate number of channels
nCh=size(d,2)/2;

% Calculate the distance between source and detector positions for each channel
for j=1:nCh
    SDdist(j,1)= sqrt((SD.SrcPos(SD.MeasList(j,1),1)-SD.DetPos(SD.MeasList(j,2),1))^2 + (SD.SrcPos(SD.MeasList(j,1),2)-SD.DetPos(SD.MeasList(j,2),2))^2 + (SD.SrcPos(SD.MeasList(j,1),3)-SD.DetPos(SD.MeasList(j,2),3))^2);
    ChPos(j,:)=(SD.SrcPos(SD.MeasList(j,1),:)+SD.DetPos(SD.MeasList(j,2),:))./2;
end

% Plot the distances between source and detector positions for each channel
%figure
%plot(SDdist)

%% Check probe design %%
SD.SrcPos=SD.SrcPos;
SD.DetPos=SD.DetPos;

% Find active measurement pairs
lst=find(SD.MeasList(:,1)>0);
ml=SD.MeasList(lst,:);
lstML = find(ml(:,4)==1);

% Plot the probe design
screenSize = get(0,'ScreenSize');
figure('Position', [ (screenSize(3)-560)/2 (screenSize(4)-560)/2 560 560])

for ii=1:nCh
    h = line( [SD.SrcPos(ml(lstML(ii),1),1) SD.DetPos(ml(lstML(ii),2),1)], ...
        [SD.SrcPos(ml(lstML(ii),1),2) SD.DetPos(ml(lstML(ii),2),2)] );
    set(h,'color',[1 1 1]*.85);
    set(h,'linewidth',4);
end
hold on
for ii=1:size(SD.SrcPos,1)
    scatter(SD.SrcPos(ii,1), SD.SrcPos(ii,2),80, 'r','filled')
    text(SD.SrcPos(ii,1), SD.SrcPos(ii,2)+0.6,['S' num2str(ii)],'color','r')
end
for ii=1:size(SD.DetPos,1)
    scatter(SD.DetPos(ii,1), SD.DetPos(ii,2),80, 'b','filled')
    text(SD.DetPos(ii,1), SD.DetPos(ii,2)+0.6,['D' num2str(ii)],'color','b')
end
for ii=1:size(ChPos,1)
    scatter(ChPos(ii,1), ChPos(ii,2),80, 'k','filled')
    text(ChPos(ii,1), ChPos(ii,2)+0.3,['Ch' num2str(ii)],'color','k')

text(ChPos(ii,1), ChPos(ii,2)-0.3,['d=' num2str(SDdist(ii))],'color','k')

end
xlim([min(min(min(SD.DetPos(:,1)),min(SD.SrcPos(:,1))),min(min(SD.DetPos(:,2)),min(SD.SrcPos(:,2))))-2 max(max(max(SD.DetPos(:,1)),max(SD.SrcPos(:,1))),max(max(SD.DetPos(:,2)),max(SD.SrcPos(:,2))))+2])
ylim([min(min(min(SD.DetPos(:,1)),min(SD.SrcPos(:,1))),min(min(SD.DetPos(:,2)),min(SD.SrcPos(:,2))))-2 max(max(max(SD.DetPos(:,1)),max(SD.SrcPos(:,1))),max(max(SD.DetPos(:,2)),max(SD.SrcPos(:,2))))+2])
set(gca, 'TickLength',[0 0])
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
title('Check Probe Design')

saveas(gca,[filename_nirs(1:end-5) '_ProbeDesign.png'])

%% RUN PREPROCESSING %%

dod = hmrIntensity2OD( d ); % Convert into OD

SD.MeasListAct=ones(Data.nCh*2,1);

%% Motion correction with wavelet
close all
iqr = 0.8
turnon=1;
[dodWavelet] = hmrMotionCorrectWavelet(dod,SD,iqr,turnon);

idx=1;
figure
for Ch=1:5
   subplot(4,2,idx)
   plot(t,dod(:,Ch),'k')
   hold on
   plot(t,dod(:,Ch+Data.nCh),'b')
   plot(t,dodWavelet(:,Ch),'m')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'r')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

idx=1;
figure
for Ch=6:9
    subplot(4,2,idx)
   plot(t,dod(:,Ch),'k')
   hold on
   plot(t,dod(:,Ch+Data.nCh),'b')
    plot(t,dodWavelet(:,Ch),'m')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'r')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

idx=1;
figure
for Ch=10:13
    subplot(4,2,idx)
   plot(t,dod(:,Ch),'k')
   hold on
   plot(t,dod(:,Ch+Data.nCh),'b')
    plot(t,dodWavelet(:,Ch),'m')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'r')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

idx=1;
figure
for Ch=14:18
    subplot(4,2,idx)
   plot(t,dod(:,Ch),'k')
   hold on
   plot(t,dod(:,Ch+Data.nCh),'b')
    plot(t,dodWavelet(:,Ch),'m')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'r')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

%% PIPELINE 1 %%
close all
% Band pass filter
fs=1/t(2);
[dodFilt,ylpf] = hmrBandpassFilt( dodWavelet, Data.fs, 0.01, 0.5 );

idx=1;
figure
for Ch=1:5
   subplot(4,2,idx)
   plot(t,dodWavelet(:,Ch),'k')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'b')
    plot(t,dodFilt(:,Ch),'r')
   hold on
   plot(t,dodFilt(:,Ch+Data.nCh),'c')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

idx=1;
figure
for Ch=6:9
    subplot(4,2,idx)
   plot(t,dodWavelet(:,Ch),'k')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'k')
    plot(t,dodFilt(:,Ch),'r')
   hold on
   plot(t,dodFilt(:,Ch+Data.nCh),'r')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

idx=1;
figure
for Ch=10:13
    subplot(4,2,idx)
   plot(t,dodWavelet(:,Ch),'k')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'k')
    plot(t,dodFilt(:,Ch),'r')
   hold on
   plot(t,dodFilt(:,Ch+Data.nCh),'r')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

idx=1;
figure
for Ch=14:18
   subplot(4,2,idx)
   plot(t,dodWavelet(:,Ch),'k')
   hold on
   plot(t,dodWavelet(:,Ch+Data.nCh),'k')
   plot(t,dodFilt(:,Ch),'r')
   hold on
   plot(t,dodFilt(:,Ch+Data.nCh),'r')
   title(['Ch ' num2str(Ch)])
   idx=idx+1;
end

% Conversion from optical density to haemoglobin conc
ppf=[5.5 4.7]; % Child

procResult.dc = hmrOD2Conc( dodFilt, SD, ppf );

%% Save info on events (in samples)

% Collaboration
idx=[];
idx=find(s(:,1)==1);
Start_Collaboration=idx(1:2:end);
End_Collaboration=idx(2:2:end);
Duration_Collaboration=End_Collaboration-Start_Collaboration;

disp(['%%%% Number of Collab Start= ' num2str(length(Start_Collaboration)) ' and End= ' num2str(length(End_Collaboration)) ' %%%%'])

% CollaborationScreen
idx=[];
idx=find(s(:,2)==1);
Start_CollaborationScreen=idx(1:2:end);
End_CollaborationScreen=idx(2:2:end);
Duration_CollaborationScreen=End_CollaborationScreen-Start_CollaborationScreen;

disp(['%%%% Number of CollabScreen Start= ' num2str(length(Start_CollaborationScreen)) ' and End= ' num2str(length(End_CollaborationScreen)) ' %%%%'])

% Individual
idx=[];
idx=find(s(:,3)==1);
Start_Individual=idx(1:2:end);
End_Individual=idx(2:2:end);
Duration_Individual=End_Individual-Start_Individual;

disp(['%%%% Number of Individual Start= ' num2str(length(Start_Individual)) ' and End= ' num2str(length(End_Individual)) ' %%%%'])

% Rest
idx=[];
idx=find(s(:,4)==1);
Start_Rest=idx(1:2:end);
End_Rest=idx(2:2:end);
Duration_Rest=End_Rest-Start_Rest;

disp(['%%%% Number of Rest Start= ' num2str(length(Start_Rest)) ' and End= ' num2str(length(End_Rest)) ' %%%%'])

%%
save([filename_nirs(1:end-5) '_ppr.mat'])

%% Plot concentrations

limMax=3*10^-6;
limMin=-3*10^-6;

% Left PFC 
idx=1;
figure
for i=1:5
    subplot(5,1,idx)
    plot(t,procResult.dc(:,1,i),'r') 
    hold on
    plot(t,procResult.dc(:,2,i),'b')
    hold on
   
    ylim([limMin limMax])
        
    xlabel('Time (s)')
    ylabel('Conc.')
 
    
    xlim([t(1) t(end)])
    title(['Ch. ' num2str(i)])
    legend('HbO2','HbR')
    idx=idx+1;
end
set(gcf,'color','w')
saveas(gca,[filename_nirs(1:end-5) '_Conc_LeftPFC_Ch1to5.png'])

% Right PFC
idx=1;
figure
for i=6:9
    subplot(4,1,idx)
    plot(t,procResult.dc(:,1,i),'r') 
    hold on
    plot(t,procResult.dc(:,2,i),'b')
    hold on
    
    ylim([limMin limMax])
        
    xlabel('Time (s)')
    ylabel('Conc.')

    xlim([t(1) t(end)])
    title(['Ch. ' num2str(i)])
    legend('HbO2','HbR')
    idx=idx+1;
end
set(gcf,'color','w')
set(gcf,'position',[106   225   735   686])
saveas(gca,[filename_nirs(1:end-5) '_Conc_RightPFC_Ch6to9.png'])

% Left TPJ
idx=1;
figure
for i=10:13
    subplot(4,1,idx)
    plot(t,procResult.dc(:,1,i),'r') 
    hold on
    plot(t,procResult.dc(:,2,i),'b')
    hold on

    ylim([limMin limMax])
        
    xlabel('Time (s)')
    ylabel('Conc.')
    
    xlim([t(1) t(end)])
    title(['Ch. ' num2str(i)])
    legend('HbO2','HbR')
    idx=idx+1;
end
set(gcf,'color','w')
set(gcf,'position',[106   225   735   686])
saveas(gca,[filename_nirs(1:end-5) '_Conc_LeftTPJ_Ch10to13.png'])

% Right TPJ
idx=1;
figure
for i=14:18
    subplot(5,1,idx)
    plot(t,procResult.dc(:,1,i),'r') 
    hold on
    plot(t,procResult.dc(:,2,i),'b')
    hold on
    
    ylim([limMin limMax])
        
    xlabel('Time (s)')
    ylabel('Conc.')
    
    xlim([t(1) t(end)])
    title(['Ch. ' num2str(i)])
    legend('HbO2','HbR')
    idx=idx+1;
end
set(gcf,'color','w')
set(gcf,'position',[106   225   735   686])
saveas(gca,[filename_nirs(1:end-5) '_Conc_RightTPJ_Ch14to18.png'])

%% (1) Visual check  of raw voltage data (Time domain and Frequency domain)

% Left PFC
figure
idx=1;
for i=1:5
    subplot(5,4,idx)
   [AX,H1,H2] =  plotyy(t,d(:,i),t,d(:,i+Data.nCh)) % Wavelength 1
    hold on
    xlabel('Time (s)')
    ylabel('Intensity (a.u.)')
    xlim(AX(1),[t(1) t(end)])
    xlim(AX(2),[t(1) t(end)])
    title(['Ch. ' num2str(i)])
    set(AX,{'ycolor'},{'k';'m'})      
    set(H1,'color','k');
    set(H2,'color','m');
    
   subplot(5,4,idx+1)
   [pxx,f] = pwelch(d(:,i),round(100*Data.fs),[],[],Data.fs); % Wavelength 1
    plot(f,20*log10(pxx),'r')
    hold on
    [pxx2,f] = pwelch(d(:,i+Data.nCh),round(100*Data.fs),[],[],Data.fs); % Wavelength 2
    plot(f,20*log10(pxx),'k')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    subplot(5,4,idx+2)
    plot(f,20*log10(pxx2),'m')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    idx=idx+4;    
end
set(gcf,'color','w')
set(gcf,'position',[1 41 1536 750])
saveas(gca,[filename_nirs(1:end-5) '_Raw_LeftPFC_Ch1to5.png'])

% Right PFC
figure
idx=1;
for i=6:9
    subplot(4,3,idx)
   [AX,H1,H2] =  plotyy(t,d(:,i),t,d(:,i+Data.nCh)) % Wavelength 1
    hold on
    xlabel('Time (s)')
    ylabel('Intensity (a.u.)')
    xlim(AX(1),[t(1) t(end)])
    xlim(AX(2),[t(1) t(end)])
    title(['Ch. ' num2str(i)])
    set(AX,{'ycolor'},{'k';'m'})      
    set(H1,'color','k');      set(H2,'color','m');
    
   subplot(4,3,idx+1)
   [pxx,f] = pwelch(d(:,i),round(100*Data.fs),[],[],Data.fs); % Wavelength 1
    plot(f,20*log10(pxx),'r')
    hold on
    [pxx2,f] = pwelch(d(:,i+Data.nCh),round(100*Data.fs),[],[],Data.fs); % Wavelength 2
    plot(f,20*log10(pxx),'k')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    subplot(4,3,idx+2)
    plot(f,20*log10(pxx2),'m')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    idx=idx+3;    
end
set(gcf,'color','w')
set(gcf,'position',[1 41 1536 750])
saveas(gca,[filename_nirs(1:end-5) '_Raw_RightPFC_Ch6to9.png'])

% Left TPJ
figure
idx=1;
for i=10:13
    subplot(6,3,idx)
   [AX,H1,H2] =  plotyy(t,d(:,i),t,d(:,i+Data.nCh)) % Wavelength 1
    hold on
    xlabel('Time (s)')
    ylabel('Intensity (a.u.)')
    xlim(AX(1),[t(1) t(end)])
    xlim(AX(2),[t(1) t(end)])
    title(['Ch. ' num2str(i)])
    set(AX,{'ycolor'},{'k';'m'})      
    set(H1,'color','k');      set(H2,'color','m');
    
   subplot(6,3,idx+1)
   [pxx,f] = pwelch(d(:,i),round(100*Data.fs),[],[],Data.fs); % Wavelength 1
    plot(f,20*log10(pxx),'r')
    hold on
    [pxx2,f] = pwelch(d(:,i+Data.nCh),round(100*Data.fs),[],[],Data.fs); % Wavelength 2
    plot(f,20*log10(pxx),'k')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    subplot(6,3,idx+2)
    plot(f,20*log10(pxx2),'m')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    idx=idx+3;    
end
set(gcf,'color','w')
set(gcf,'position',[1 41 1536 750])
saveas(gca,[filename_nirs(1:end-5) '_Raw_LeftTPJ_Ch10to13.png'])

% Right TPJ
figure
idx=1;
for i=14:18
    subplot(5,4,idx)
   [AX,H1,H2] =  plotyy(t,d(:,i),t,d(:,i+Data.nCh)) % Wavelength 1
    hold on
    xlabel('Time (s)')
    ylabel('Intensity (a.u.)')
    xlim(AX(1),[t(1) t(end)])
    xlim(AX(2),[t(1) t(end)])
    title(['Ch. ' num2str(i)])
    set(AX,{'ycolor'},{'k';'m'})      
    set(H1,'color','k');      set(H2,'color','m');
    
   subplot(5,4,idx+1)
   [pxx,f] = pwelch(d(:,i),round(100*Data.fs),[],[],Data.fs); % Wavelength 1
    plot(f,20*log10(pxx),'r')
    hold on
    [pxx2,f] = pwelch(d(:,i+Data.nCh),round(100*Data.fs),[],[],Data.fs); % Wavelength 2
    plot(f,20*log10(pxx),'k')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    subplot(5,4,idx+2)
    plot(f,20*log10(pxx2),'m')
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
    xlim([f(1) 6])
    title(['Ch. ' num2str(i)])
    
    idx=idx+4;    
end
set(gcf,'color','w')
set(gcf,'position',[1 41 1536 750])
saveas(gca,[filename_nirs(1:end-5) '_Raw_RightTPJ_Ch14to18.png'])

