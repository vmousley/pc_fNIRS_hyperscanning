%% Config file for WTC Pipeline
%{

Written by Victoria Mousley, Centre for Brain and Cognitive Development
Email: v.mousley@bbk.ac.uk

Set your study-specific paths & variables to run the WTC pipeline. 

%}

% Paths

% analysisPath = '\set\your\path'; % <<<<<<< SET
% basePath = '\set\your\path\'; % <<<<<<< SET
% exportPAth = '\set\your\path'; % <<<<<<< SET

% Number of channels
chanNum = 18;

% Channels belonging to each region of interest
leftPFC = [1 2 4 5];
rightPFC = [6 7 8 9];
leftTPJ = [10 11 12 13];
rightTPJ = [14 15 17 18];

% Channel numbers corresponding to short-separation channels
sscs = [3; 16];

% Frequency band of interest
poi = [10 50]; 

% Analysis window for expeirmental conditions (seconds)
cond_len = 120;

% Analysis window for rest conditions (seconds)
rest_len = 80; 

% NIRS sampling rate
samp = 25;


