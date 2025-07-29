
function DataFinal = traces(Data)


DataFinal = Data;

FileName = Data.FileName;
FilePath = Data.FilePath;

MaxDistancebig = Data.Traces.MaxDistancebig; % Maximum distance (Pixels) between objects in consecutive frames
MinLength = Data.Traces.MinLength;
FigNumber = 2;

Position = Data.Position;

tic
XB = Position(1).X1;   %Coordinates in the 1st frame
YB = Position(1).Y1;


% Create a set of trace structures using the initial cooordinates and setting t=0 
 for j = 1:1:length(XB)              
    Tracebig(j).TB = [0];
    Tracebig(j).XB = [XB(j)];
    Tracebig(j).YB = [YB(j)];
 end

%beginframe=50;
for i = 2:1:length(Position)
    disp(['** TRACES (' FileName ')- position ' int2str(i) '/' int2str(length(Position)) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])    
    XB = Position(i).X1; %Coordinates in frame i (e.g, frame 2)
    YB = Position(i).Y1;
      
%% Big particles
    
    for p = 1:1:length(Tracebig)
        Distancebig = sqrt((XB-Tracebig(p).XB(end)).^2 + (YB-Tracebig(p).YB(end)).^2 ); %Check distance 
        MinDistanceIndexbig = find(Distancebig==min(Distancebig));
        if (length(MinDistanceIndexbig)>0)
            MinDistanceIndexbig = MinDistanceIndexbig(1);
            if (Distancebig(MinDistanceIndexbig)<MaxDistancebig) && (Tracebig(p).TB(end) == (i-2))    
                Tracebig(p).TB = [Tracebig(p).TB (i-1)];
                Tracebig(p).XB = [Tracebig(p).XB XB(MinDistanceIndexbig)];
                Tracebig(p).YB = [Tracebig(p).YB YB(MinDistanceIndexbig)];
                XB(MinDistanceIndexbig) = Inf;
                YB(MinDistanceIndexbig) = Inf;
            end
        end
    end
    
    for k = 1:1:length(XB)
        if (XB(k)<Inf && YB(k)<Inf) 
            p = p+1;
            Tracebig(p).TB = [(i-1)];
            Tracebig(p).XB = [XB(k)];
            Tracebig(p).YB = [YB(k)];
        end
    end
end


i=0;


DataFinal.Traces.Tracebig = Tracebig;


 figure(1)
    set(gcf,'Units','normalized','Position',[0 0 1 1])
    axes('Position',[.05 .1 .8 .8])
    hold on
    for j = 1:1:length(Tracebig)
        plot(Tracebig(j).XB,Tracebig(j).YB,'r')
    end
    hold on
  

    hold off
    axis equal
    xlabel('Pixels')
    ylabel('Pixels')
    axes('Position',[.88 .1 .1 .8])
    hold on
    for j = 1:1:length(Tracebig)
        plot(Tracebig(j).TB,j*ones(size(Tracebig(j).TB)),'r')
    end
    hold on
    hold off
    xlabel('Frames')

    colormap(bone)
    

end
