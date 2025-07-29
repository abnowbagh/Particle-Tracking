clear all;
close all;
clc ;


[DATA.File, DATA.Path] = uigetfile('*.mat');
load([DATA.Path DATA.File]);

DATA.Trajectories1 = Data.Trajectories.Trajectorybig;
DATA.Position = Data.Position; 

cal=0.293; %um/pixel
fps=5;

%% MSD
%Enter calibration factor (pixel to micron) and frame rate
% cal=1 if X,Y are already in physical units

msd_mean=NaN(length(DATA.Trajectories1),length(DATA.Position)-1);
msd_std=NaN(length(DATA.Trajectories1),length(DATA.Position)-1);
v_mean=NaN(length(DATA.Trajectories1),length(DATA.Position)-1);
v_std=NaN(length(DATA.Trajectories1),length(DATA.Position)-1);
time=NaN(length(DATA.Trajectories1),length(DATA.Position)-1);

v_mean_par=[]; %for mean velocity per particle (if you don't want to fit)
msd_mean_par=[];

for i = 1:1:length(DATA.Trajectories1) %going through the list of trajectories
    TB = DATA.Trajectories1(1,i).TB;
    XB = DATA.Trajectories1(1,i).XB;
    YB = DATA.Trajectories1(1,i).YB;        
    tic
    diff1 = NaN(max(TB)-1,max(TB)-1); % define a "place-holder"
    diff2 = NaN(max(TB)-1,max(TB)-1); % define a "place-holder"
    msd1 = 0;
     %counter for time

    for k = 2:length(TB)-1 % "from"
        for l = 1:k-1 % "to"
            % define time, this has to be done in case frames are skipped
            T1 = TB(k)-TB(1); % "from"
            T2 = TB(l)-TB(1); % "to"
            msd1 = (XB(k)-XB(l))^2+(YB(k)-YB(l))^2;
            disp1 = sqrt((XB(k)-XB(l))^2+(YB(k)-YB(l))^2);
            diff1(T1,T1-T2) = msd1; %MSD for every possible delta T (frame vs delta T)
            diff2(T1,T1-T2) = disp1;
        end
    end
    diff1 = diff1*cal.^2;
    diff2 = diff2*cal;
    delta=NaN(max(TB)-1,1)';
    delta(1:max(TB)-1)=[1:1:max(TB)-1];
    msd_mean(i,1:(max(TB)-1)) = nanmean(diff1);
    msd_std(i,1:(max(TB)-1))  = nanstd(diff1);
    v_mean(i,1:(max(TB)-1)) = nanmean(diff2)./delta;
    v_mean(i,1:(max(TB)-1)) = v_mean(i,1:(max(TB)-1))*fps;
    v_std(i,1:(max(TB)-1))  = nanstd(diff2)*fps./delta;
    
    if (nanmean(v_mean(i,:))>0.1)
        v_mean_par= [nanmean(v_mean(i,:)), v_mean_par];
        msd_mean_par=[nanmean(msd_mean(i,:)), msd_mean_par];
    end
    
    %time=TB(1:length(msd_mean(i,:)))./5;
    msd_nan=msd_mean(i,:);
    msd_clean=msd_nan(~isnan(msd_nan));
    %time=TB(2:end-1)./5;

    figure(i)
    subplot(2,1,1)
    plot(msd_mean(i,:));
    title('MSD')
    xlabel('frame [5 fps]')
    ylabel('MSD')
    
    v_nan=v_mean(i,:);
    v_clean=v_nan(~isnan(v_nan));    
    subplot(2,1,2)
    plot(v_mean(i,:));
    title('Velocity')
    xlabel('frames [1 fps]')
    ylabel('v')
    
    %saveas(gcf,['msd and velocity' i '.jpg'])
    time=[];
    
    
end
% avg= mean(v_mean_par);
% std= std(v_mean_par);
% disp('average velocity and std');
% disp (avg);
% disp(std);

% MSD=nanmean(msd_mean);
% DELTA_TIME=[1:1:length(MSD)]./fps;
% figure(1)
% plot(DELTA_TIME,msd_mean)
% xlabel('Time (s)')
% ylabel('MSD (um^2)')

%%  fitting with D constant
% % Curve Fitting Parameters 1(only v )
% 
% 
% t_max = 49; % t_max < tau_r
% t = linspace(0, t_max, t_max + 1);
% D = 0.086; % Initial guess for Diffusion coefficient
% V = 1; % Initial guess for Velocity
% vel=[]; %to store the velocities
% 
% % Modified fitting function with only V as a parameter
% func = @(p, t) 4 * D * t + p^2 * t.^2;
% 
% 
% for i=1:1:size(msd_mean,1)
%     m_i=msd_mean(i,:); %getting first t_max msd values
%     msd_i= [0, m_i(1:t_max)]; %corresponding to t=0
%     hasNan= any(isnan(msd_i));
%     
%     % % Calculate the average MSD for all trajectories
%     % avg_msd = mean(msd_mean(:, 1:t_max), 1);
%     % avg_msd = [0, avg_msd]; % Adding 0 for t=0
%     
%     % Curve Fitting for MSD data
%     if(~hasNan)
%     initialguess = V; % Only one parameter now
%     lb = []; % Lower bound for V (can be empty)
%     ub = []; % Upper bound for V (can be empty)
%     [p_fit, resnorm, residual, exitflag, output] = lsqcurvefit(func, initialguess, t, msd_i, lb, ub);
%     
%     % Display the fitted parameter V (D is fixed)
%     %disp('Fitted parameter V:');
%     %disp(p_fit);
%     
%     %Calculating R_adj
%     y_fit = 4 * D * t + p_fit(1)^2 * t.^2;
%     SStot= sum((msd_i-mean(msd_i)).^2); % SSE
%     SSe= sum((msd_i-y_fit).^2); %SST
%     adj=(length(msd_i)-1)/(length(msd_i)-3); %adjusting factor for R
%     Radj=1-adj*SSe/SStot; %R^2 adjusted
%     
%     if Radj>0.96
%         vel=[p_fit, vel];
%         % Plot the original data and the fitted curve
%         figure;
%         scatter(t, msd_i, 'o', 'DisplayName', ' MSD Data');
%         hold on;
%         plot(t, y_fit, 'r-', 'DisplayName', 'Fitted Curve');
%         xlabel('t');
%         ylabel('msd');
%         legend('show');
%         title('Fitting y = 4Dt + V^2t^2 for MSD');
%     end
%     end
% end
% disp('Average Velocity and Stdev');
% disp(mean(vel));
% disp(std(vel));
% %

%%
% %fitting using polyfityero (both D and v)
% v0=[];
% D_t=[];
% t_max=80; %t_max< tau_r
% t =linspace(0,t_max,t_max+1);
% t1=linspace(0,t_max);
% coeff=zeros(4,size(msd_mean,1));
% figure
% 
% D=0.086;
% func=@(v,t) 4*D.*t*v^2*t.^2;
% 
% for i=1:1:size(msd_mean,1)
%     m_i=msd_mean(i,:); %getting first t_max msd values
%     msd_i= [0, m_i(1:t_max)]; %corresponding to t=0
%     
%     initialguess = [1];
%     [v_fit, resnorm, residual, exitflag, output] = lsqcurvefit(func, initialguess, msd_i, t);
%     
%     Display the fitted parameters
%     disp('Fitted parameter v:');
%     disp(v_fit);
%     
%     Plot the original data and the fitted curve
%     figure;
%     scatter(t, msd_i, 'o', 'DisplayName', 'Data');
%     hold on;
%     
%     y_fit = 4*D*t1+v_fit^2*t1.^2;
%     plot(t1, y_fit, 'r-', 'DisplayName', 'Fitted Curve');
%     xlabel('t');
%     ylabel('msd');
%     legend('show');
%     title('Fitting y = 4Dt + v^2x^2');
% end
% % %     
% %     
% %     [p,S]=polyfitZero(t, msd_i,2); %fitting with 0 intercept
% %     
% %     msd_fit=polyval(p,t); %fit at discrete points
% %     msd_curve=polyval(p,t1); %fit at continuous (more) points
% %     SStot= sum((msd_i-mean(msd_i)).^2); % SSE
% %     SSe= sum((msd_i-msd_fit).^2); %SST
% %     adj=(length(msd_i)-1)/(length(msd_i)-3); %adjusting factor for R
% %     Radj=1-adj*SSe/SStot; %R^2 adjusted
% %     Rsq=1-SSe/SStot; %R^2
% %     
% %     if(p(1))>0 && Radj>0.9 %selecting only values which fit perfectly
% %         vel=sqrt(p(1));
% %         coeff(1,i)=vel;
% %         D=p(2)/4; 
% %         coeff(2,i)=D;
% %         coeff(3,i)=Rsq;
% %         coeff(4,i)=Radj;
% %         v0=[vel v0];
% %         D_t=[D D_t];
% %             figure
% %     plot(t,msd_i,'o')
% %     hold on
% %     plot(t1,msd_curve)
% %    
% %     txt =['v: ' num2str(vel) ' D: ' num2str(D) ' R2: ' num2str(Rsq)];
% %     text (10,3000,txt)
% %     hold off
% %     
% %     end
% %    
% % 
% %         
% %    
% %     
% %     
% %     
% %     
% %     figure
% %     plot(t,msd_i,'o')
% %     hold on
% %     plot(t1,msd_curve)
% %    
% % %     txt =['v: ' num2str(vel) ' D: ' num2str(D) ' R2: ' num2str(Rsq)];
% % %     text (10,3000,txt)
% %     hold off
% %    
% %    
% % end
% % 
% % 
% % velociy=mean(v0)
% % std_v0=std(v0)
% % Diffusion = mean(D_t)
% % std_Diff= std(D_t)
% % % t1_2=t1.^2;
% % % y=4*mean(D_t)*t1+(mean(v0)^2)*t1_2;
% % % %plot(t1,y)

%%

    
    
    
    




