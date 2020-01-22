%% script to plot slip history stairs from a .out file
%
% Requires sliphistory.out/slips file 
% User must input the site specific details for the SL code in the section
% below
%
%
% LCG Nov 2018

clear
close all
%% Set up for plotting fit to data

faultname = 'Caprociano';

slip_file = '../sliphistory.out';

burnin = 40000;

% Number to plot in parula
top_slips = 20;

grey_slips = 100;


%% No further input required below this point

% Load sliphistory file
SlipHistories = dlmread(slip_file);

% Read offsets
slip = SlipHistories(1,4:end);

% Exclude top line and burnin from file
SlipHistories = SlipHistories(burnin:end,:);
SizeOfFile = size(SlipHistories)
maxmodel = SizeOfFile;

% Maxlik calculation
% Sort Matrix based on likelihood
    [values, order] = sort(SlipHistories(:,1));
    SortedSlipHistories = SlipHistories(order,:);
    SortedSlipHistories = SortedSlipHistories(:,4:end);
    Top1000 = SortedSlipHistories(end-top_slips:end,:);
    Top9000 = SortedSlipHistories(end-grey_slips:end,:);
    

%% Slip vector

NumberOfModels = size(SortedSlipHistories);

 % Calculate cumulative slip vector
     CumSlip = slip;
     CumSlip(length(CumSlip)) = -slip(length(slip));
     CumSlip(length(CumSlip)-1) = 0;
     for i = (length(CumSlip)-2):-1:1
         CumSlip(i) = slip(i) + CumSlip(i+1);
     end
     
     % Create vector of zeros of same length as slip vector
     CumSlip = zeros(1,length(slip));
 
     % First entry should be - trench depth
     CumSlip(length(CumSlip)) = -slip(length(slip));
 

     % Add all slip events
     for i = (length(CumSlip)):-1:2
         CumSlip(i-1) = slip(i) + CumSlip(i);
     end
     
     %Add one more event at scarp age
     CumSlip	= [(CumSlip(1)+slip(1)),CumSlip];
     
%% Ploting

    % Figure
    figure(32);
    hold on;
    
    ax = gca;
    ax.XRuler.Exponent = 0;
    
    xlabel('Time (yrs)');
    ylabel('slip (cm)');
    title(sprintf('Best fit slip histories: %s', faultname));
    
    %Set specific axes here (e.g. for all plots to be the same)
    ylim([0 800])
    xlim([0 25000])
    
    % Plot slip histories in color scale according to order in file
    % (likelihood) - flipped so the 'best fits' are on top
    
%      cc = parula(NumberOfModels(1)); 
%      cc = flipud(cc);
%      SortedSlipHistories = flipud(SortedSlipHistories);
%      
%      for i = 1:NumberOfModels(1)
%          ages = SortedSlipHistories(i,:);
%          ages = [ages(1),ages];
%          stairs(ages,CumSlip,'color',cc(i,:));
%      end

% Plot slip histories with top 1000 in parula, rest in grey

    %agestoplot = SortedSlipHistories(1:grey_slips,:);
    
    for k = 1:grey_slips
        agesgrey = Top9000(k,:);
        agesgrey = [agesgrey(1),agesgrey];
        stairs(agesgrey,CumSlip,'Color',[0.5 0.5 0.5]);
    end
    
    hold on;

      cc = parula(top_slips); 
      cc = flipud(cc);
     Top1000 = flipud(Top1000);
     
     for i = 1:top_slips
         ages = Top1000(i,:);
         ages = [ages(1),ages];
         stairs(ages,CumSlip,'color',cc(i,:));
     end
     
 
     
%% Plot most likely slip history on top - not currently coded in but could 
%  be added

%age_plot = [age(1),age];
%stairs(age_plot,CumSlip,'red','LineWidth',1.8);
 
% ylim([0 maxheight]) %new for plotting all at same size
% xlim([0 23000])

% Add another slip history plotted on top of histoplot manually, using same
% slips as in maxlik file: 

% ageC = [15500 13562 11625 9687 7750 5812 3875 1937 0];
% ageC = [ageC(1),ageC];
% stairs(ageC,CumSlipM,'blue','LineWidth',1.8);
