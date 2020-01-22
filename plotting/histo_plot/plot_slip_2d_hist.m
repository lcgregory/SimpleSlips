%% script to plot 2D histograms of slip histories 
%
% This script makes a nice plot of the density of your slip histories. It
% does not necessarily reflect the detailed pattern of slip, because two
% different histories may pass through the same point on the plot. Compare
% this plot with the top 10k models to see the detail of your fits. This
% plot is good for showing the statistical distribution of slip (and will
% correspond to the slip rate plots) but it is not good for seeing the
% tradeoff between scarp age and slip rate.
%
% If more examples of the best fit/constant rateslip histories are to be 
% plotted, they must be copied in the plotting section. 
%
% TJW and LC 6 Feb 2018. Modified Nov 2018.
clear
close all

%% Enter number of slips in sliphistory.out and a burnin on the model
%maxheight = 820; %Does not include trench
burnin = 40000;
tn = 200; %number of time bins that the history is divided into, usually 200
cd = 15000 ; %Colour density, needs to scale with number of models by ~10%
slipfile = '../../sliphistory.out';
faultname = 'Caporciano';

%% No further input required here
% Load models and sort one version
sliphistory = load(slipfile);

% Read offsets (maxlik is the same as later slip)
slipM = sliphistory(1,4:end);
% Exclude offsets and rename for sorting (keep sliphistory intact)
sorting = sliphistory(2:end,:);
% Remove burnin
sorting = sorting(burnin:end,:);
% Sort maxlik based on likelihood
[values, order] = sort(sorting(:,1));
maxlik = sorting(order,:);

%Find the size (sn) of the slips and subtract the trench
[numRows,sn] = size(slipM);
sn = sn - 1;

%Calculate the total height of the scarp above the trench (for making the
%plot at the correct height)
maxheight = sum(slipM(1,1:end-1));

%% Maxlik calculation

%take most likely age from maxlik variable
ageM = maxlik(end,4:end);

% Calculate cumulative slip vector
    CumSlipM = slipM;
    CumSlipM(length(CumSlipM)) = -slipM(length(slipM));
    CumSlipM(length(CumSlipM)-1) = 0;
    for i = (length(CumSlipM)-2):-1:1
        CumSlipM(i) = slipM(i) + CumSlipM(i+1);
    end
    
    % Create vector of zeros of same length as slip vector
    CumSlipM = zeros(1,length(slipM));

    % First entry should be - trench depth
    CumSlipM(length(CumSlipM)) = -slipM(length(slipM));

    % Add all slip events
    for i = (length(CumSlipM)):-1:2
        CumSlipM(i-1) = slipM(i) + CumSlipM(i);
    end
    
    %Add one more event at scarp age
    CumSlipM	= [(CumSlipM(1)+slipM(1)),CumSlipM];
    
    
    
    

%% For histoplot

sliphistory(:,1:3)=[];
sliphistory(:,end)=[];
slipsizes=sliphistory(1,:);

sliptimes=sliphistory(burnin:end,:);

cumslip = cumsum(slipsizes);
cumslip = cumslip(end)-cumslip;

A = cumslip(1,end-1)/2;
slips = cumslip+A;
slips = fliplr(slips);

maxheightA = cumslip(1,1);

cumslipmat=repmat(cumslip,size(sliptimes,1),1);

timescol=reshape(sliptimes,numel(sliptimes),1);
slipcol=reshape(cumslipmat,numel(sliptimes),1);
times=linspace(0,max(timescol),tn);
histfile=hist2d([timescol slipcol], tn, sn,[0 max(timescol)],[0 max(slipcol)]);

%% For synthetic data plot without the same slips

% slipsyn = [100 100 40 40 100 20 40 100 250 80];
% agesyn = [12000 5400 2208 2072 1573 1532 1398 1215 785 0];
% 
% % Calculate cumulative slip vector
%     CumSlipS = slipsyn;
%     CumSlipS(length(CumSlipS)) = -slipsyn(length(slipsyn));
%     CumSlipS(length(CumSlipS)-1) = 0;
%     for i = (length(CumSlipS)-2):-1:1
%         CumSlipS(i) = slipsyn(i) + CumSlipS(i+1);
%     end
%     
%     % Create vector of zeros of same length as slip vector
%     CumSlipS = zeros(1,length(slipsyn));
% 
%     % First entry should be - trench depth
%     CumSlipS(length(CumSlipS)) = -slipsyn(length(slipsyn));
% 
%     % Add all slip events
%     for i = (length(CumSlipS)):-1:2
%         CumSlipS(i-1) = slipsyn(i) + CumSlipS(i);
%     end
%     
%     %Add one more event at scarp age
%     CumSlipS	= [(CumSlipS(1)+slipsyn(1)),CumSlipS];
    
%% For plotting

%subplot(1,2,1)
imagesc(times,slips,histfile)

ax = gca;
ax.XRuler.Exponent = 0;

title(sprintf('Density of all modelled slip histories: %s',faultname));
ylim([0 maxheight]) %plotting all at same size
xlim([0 24000]) 
xlabel('Time (yrs)');
ylabel('slip (cm)');

axis xy
colormap(flipud(parula))
c=colormap;
c(1,:)=[1 1 1];
colormap(c)
h = colorbar('vert', 'east');
ylabel(h,'Number of models');

caxis([0 cd])

hold on
% Plot highest likelihood slip history on top of histoplot in red
ageM = [ageM(1),ageM];
stairs(ageM,CumSlipM,'r','LineWidth',1.2);
%stairs(ageM,CumSlipM,'Color',[0.5 0.5 0.5],'LineWidth',1.8);

% Add another slip history plotted on top of histoplot manually, using same
% slips as in maxlik file: 

%    agesyn = [agesyn(1),agesyn];
%    stairs(agesyn,CumSlipS,'red','LineWidth',1.8);


%% cumulative sliprate figure

% calculate cumulative slip rate - this is averaged with time over the
% duration of slip in the past - e.g. at 10 kyr the rate is fully averaged
% over 10 kyr based on the height at that time

figure(2)

slipcolmm = slipcol*10;
ratecol = slipcolmm./timescol;

%rn is the 'rate number', which corresponds to the number of divisions for
%the slip rate (100 is standard). 
%tn=200;
timesrate=linspace(0,max(timescol),tn);
rn=100;
rate=linspace(0,max(ratecol),rn);
histfileR=hist2d([timescol ratecol], tn, rn,[0 max(timescol)],[0 max(ratecol)]);
%subplot(1,2,1)
imagesc(timesrate,rate,histfileR)


%ylim([0 800]) 
xlabel('Time (yrs)');
%ylabel('Slip (cm)');
title(sprintf('Time averaged slip rate (mm/yr): %s',faultname));

axis xy
colormap(flipud(parula))
c=colormap;
c(1,:)=[1 1 1];
colormap(c)
colorbar('vert', 'east')
ylim([0 4])
caxis([0 cd])

