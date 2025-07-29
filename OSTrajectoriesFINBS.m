clear all;
close all; 
clc;

%% Chose the file

[DATA.File,DATA.Path, DATA.Frame] = uigetfile('*.mat');

disp(['** File: ' DATA.File])
disp(['** Path: ' DATA.Path])

Data.FileName = DATA.File;
Data.FilePath = DATA.Path;
load([DATA.Path DATA.File]);
Data.Position = DATA.Frame; 


% TRACES
FileTraces = cat(2,Data.FilePath,Data.FileName(1:end-4),'_traces');

if exist([FileTraces '.mat']) 
    load([FileTraces '.mat'])
    disp(['** Loaded Traces - ' FileTraces])
else
    Data.Traces.MaxDistancebig = 50; % Maximum distance (Pixels) between objects in consecutive frames (originally 20)
    Data.Traces.FigNumber = 2;
    Data.Traces.MinLength = 2;   %Minimum length of the traces
    Data = OStracesBS(Data);
%      save([FileTraces '.mat'],'Data')
%      saveas(gcf,[FileTraces '.jpg'],'jpg')
end

%% TRAJECTORIES

FileTrajectories = cat(2,Data.FilePath,Data.FileName(1:end-4),'_trajectories');
% if exist([FileTrajectories '.mat']) 
%     load([FileTrajectories '.mat'])
%     disp(['** Loaded Trajectories - ' FileTrajectories])
% else
    Data.Trajectories.MaxDistance = 30; % Maximum distance (pixels) between objects in consecutive frames (originally 5)
    Data.Trajectories.MaxTime = 5; % Maximum time that can be skipped
    Data.Trajectories.FigNumber = 2;
    Data.Trajectories.MinLength = 20; %originally 100
    Data = OStrajectoriesBS(Data);
     save([FileTrajectories '.mat'],'Data')
     saveas(gcf,[FileTrajectories '.jpg'],'jpg')
% end

