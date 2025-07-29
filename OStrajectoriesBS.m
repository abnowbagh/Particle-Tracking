function DataFinal = trajectories(Data)
DataFinal = Data;
FileName = Data.FileName;
FilePath = Data.FilePath;

MaxDistance = Data.Trajectories.MaxDistance; % Maximum distance (Pixels) between objects in consecutive frames
MaxTime = Data.Trajectories.MaxTime; % Maximum number of frames that can be skipped
MinLength = Data.Trajectories.MinLength;
FigNumber = Data.Trajectories.FigNumber;

Tracebig = Data.Traces.Tracebig;


tic
i = 0;

%% Big Particles

for k = 1:1:length(Tracebig)
    if length(Tracebig(k).XB)>MinLength
    disp(['** TRAJECTORIES (' FileName ')- trace ' int2str(i) '/' int2str(length(Tracebig)) ' - ' int2str(toc) '.' int2str(mod(toc,1)*10) 's'])    
    if (Tracebig(k).TB~=Inf)
        i = i+1;
        Trajectorybig(i).TB = Tracebig(k).TB;
        Trajectorybig(i).XB = Tracebig(k).XB;
        Trajectorybig(i).YB = Tracebig(k).YB;
        for j = 2:1:length(Tracebig)
            TEnd = Trajectorybig(i).TB(end);
            XEnd = Trajectorybig(i).XB(end);
            YEnd = Trajectorybig(i).YB(end);

            TStart = Tracebig(j).TB(1);
            XStart = Tracebig(j).XB(1);
            YStart = Tracebig(j).YB(1);

            if ( (TStart>TEnd) && ...
                    (TStart<(TEnd+MaxTime)) && ...
                    ( (((XStart-XEnd))^2+((YStart-YEnd))^2) < MaxDistance^2 ) )
                Trajectorybig(i).TB = [Trajectorybig(i).TB Tracebig(j).TB];
                Trajectorybig(i).XB = [Trajectorybig(i).XB Tracebig(j).XB];
                Trajectorybig(i).YB = [Trajectorybig(i).YB Tracebig(j).YB];
                Tracebig(j).TB = Inf;
                Tracebig(j).XB = Inf;
                Tracebig(j).YB = Inf;
            end
%         end
            end
        end
    end
end

tic
i = 0;

DataFinal.Trajectories.Trajectorybig = Trajectorybig;


if (FigNumber>0)
    color_list=['b','g','r','c','m','k'];
    color_rand=randi(6,1,length(Trajectorybig));
    figure(FigNumber)   
    set(gcf,'Units','normalized','Position',[0 0 1 1])
    axes('Position',[.05 .1 .8 .8])    
    hold on
    for j = 1:1:length(Trajectorybig)        
         plot(Trajectorybig(j).XB,Trajectorybig(j).YB,color_list(color_rand(j)))
%          plot(Trajectorybig(j).XB,Trajectorybig(j).YB,'r'))
    end
    hold on
    hold off
    axis equal
    xlabel('Pixels')
    ylabel('Pixels')
    axes('Position',[.88 .1 .1 .8])
    hold on
    for j = 1:1:length(Trajectorybig)
        plot(Trajectorybig(j).TB,j*ones(size(Trajectorybig(j).TB)),color_list(color_rand(j)))        
    end
    hold on
    hold off
	xlabel('Frames')
    colormap(bone)
end