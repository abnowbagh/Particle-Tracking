clear all;
close all; 
clc;

%% Chose the file

[Data.File,Data.Path] = uigetfile('*.*');

disp(['** File: ' Data.File])
disp(['** Path: ' Data.Path])

Video = VideoReader([Data.Path Data.File]);

D=5; %diameter of particle in microns
cal=0.293; %calibration in um/pixel (0.293 for 20x, 0.586 for 10x, 0.147 for 40x)
%% Loop on video frames

begframe=1;
endframe=Video.NumberOfFrames;
AF=[]; %array for area fractions of all ROIs

%% if video

for a =begframe:1:endframe
    Data.Frame.Number = a;
     a
     A = read(Video,a); % Original image 

% %% if image
%      A=imread ([Data.Path Data.File]);
%      %image(A)
%     %imflatfield might be good to flatten image
    
    %%
     F1 = fspecial('gaussian',10,5); % filter type, size, sigma
     B1 = imfilter(A,F1,'replicate'); % Background    
     A1 = A-B1; % Image - Background
     S1 = imsharpen(A, 'Radius',10,'Amount',10);%imadjust(A1, stretchlim(A1),[]); % Image with stretched colors

     channel1 = 1;
     level1 = .36;
    BW1 = im2bw(S1(:,:,1),level1); % Image in BW

%    BW1=imbinarize(S1(:,:,1),'global');
%     figure, imshow(BW1)
%  
%     dilate1 = strel('disk',1);
%     BW1 = imdilate(BW1,dilate1);
%     erode1 = strel('disk',1);
%     BW1 = imerode(BW1,erode1);
% %         
%      BW2 = BW1;
%     BW1 = 1 - BW1;
%     BW1 = imfill(BW1,'holes');

 %   imshow(BW1)
%     imshow(BW2)

    regions1 = bwlabel(BW1,8);
    props1 = regionprops(regions1, 'Centroid', 'Area','Eccentricity');
% %
%     regions2 = bwlabel(BW2,8);
%     props2 = regionprops(regions2, 'Centroid', 'Area','Eccentricity');

    minarea1 = 500;          % minimum size for detected area
    maxarea1 = 2000;       % minimun size for detected area

    
    X1 = []; %Coordinates of the particle
    Y1 = [];


% Particles darker than background    
% % 
    for n = 1:1:size(props1)
        if (props1(n).Area > minarea1) && (props1(n).Area < maxarea1)
            X1 = [X1; props1(n).Centroid(1)];
            Y1 = [Y1; props1(n).Centroid(2)];            
        end
    end
%     
% % % Particles brighter than background
% % % 
%       for n = 1:1:size(props2)
%           if (props2(n).Area > minarea1) && (props2(n).Area < maxarea1) %&&(props2(n).Eccentricity < 0.9)         
%             X1 = [X1; props2(n).Centroid(1)];
%             Y1 = [Y1; props2(n).Centroid(2)];
%          end       
%        end  
% % %     

%% Output

    Data.Frame.X1 = X1;
    Data.Frame.Y1 = Y1;
    Data.Frame.NumberOfParticles=length(X1);
     
      %% Plots
    %if (mod(a,200)==0)
      figure(1)
    image(A)
    hold on
    scatter(X1,Y1,16,'+b')
    %title(['Area Fraction' ])    
    title(['Frame number ' int2str(a)])  
    drawnow();
    %end
    
%     Area_particles= Data.Frame.NumberOfParticles * 3.14* (D/2)^2; %in micron sq
%     Area_image= length(A(:,1,1))*length(A(1,:,1))*cal^2; %in micron sq
%     Area_frac= Area_particles/Area_image*100;
%     disp(Area_frac)
%     AF=[AF Area_frac];
    
%    SAVE (actual)
   d = [Data.Path Data.File(1:end-4)];
     if (isdir(d)==0)
         mkdir(d)
     end
     save([Data.Path Data.File(1:end-4) '\' Data.File(1:end-4) '_psi6_' int2str(a) '.mat'],'Data')


end
% mean(AF)
% std(AF)
             
 
% Save


DATA.File = Data.File;
DATA.Path = Data.Path;
DATA.FrameRate = Video.FrameRate;

%Actual
for a = begframe:1:endframe
    load([Data.Path Data.File(1:end-4) '\' Data.File(1:end-4) '_psi6_' int2str(a) '.mat'])
    a
    DATA.Frame(a).X1 = Data.Frame.X1;
    DATA.Frame(a).Y1 = Data.Frame.Y1;
    DATA.Frame(a).NumberOfParticles = Data.Frame.NumberOfParticles;

  end
save([Data.Path Data.File(1:end-4) '.mat'],'DATA')



