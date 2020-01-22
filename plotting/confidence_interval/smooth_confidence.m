%% script to determine and plot confidence intervals 
% The confidence intervals are calculated as percentiles of the model
% distribution - e.g. showing the area covering 95% of the models. 
% This calculation is simply based on the distribution, not the fit to the
% data or liklihood function (though the model results are related to these
% two parameters, so they are linked).
%
%   HG 2018 & LG 2020
%

clear
%close all
%% Set up 

faultname = 'Caporciano';

% Direct to sliphistory.out file
slipfile = '../../sliphistory.out';

% Set burnin
burnin = 40000;

% set percentage bounds - to calculate the 95% percentile use 2.5 and 97.5
low_bound = 2.5;
high_bound = 97.5;

% set plotting options, e.g. ['k','.'] for black dots
PlotOpts_95 = 'k';
PlotOpts_avg = 'r';
PlotOpts_maxlik = 'b';

%% No further input required below this point

%Load sliphistories
SlipHistories = dlmread(slipfile);
%Extract slip magnitude
slip = SlipHistories(1,4:end-1);

slip = [slip 0];

%Extract slip ages minus the burnin and the trench
ages = SlipHistories(burnin:end,4:end-1);
ages = [ages ages(:,end)];
SizeOfFile = size(SlipHistories);

%Flip timing of slip (ts) and magnitude of slip (ms) left ro right 
ages=fliplr(ages);
slip=fliplr(slip);
    
% Calculate the 95th percentile bounds (change the [2.5] and [97.5] 
% if you want to calculate a different range
low = prctile(ages,low_bound,1);
high = prctile(ages,high_bound,1);

% Calculate the mean of the dataset
avg = mean(ages,1);

% calculate cumulative slip
slip = cumsum(slip);


% Calculate the most likely model in the distribution to plot on top
% Sort Matrix based on likelihood or RMSw
    [values, order] = sort(SlipHistories(:,1));
    SortedSlipHistories = SlipHistories(order,:);
    
% Set age array from sorted slip histories

maxlik = SizeOfFile(1) - 1;
maxlikage = SortedSlipHistories(maxlik,4:end-1);
maxlikage = [maxlikage maxlikage(:,end)];
maxlikage = fliplr(maxlikage);


%% Plot the results

% Figure
figure(9)
hold on;

ax = gca;

xlabel('Time (yrs)');
ylabel('slip (cm)');
title(sprintf('Model distribution: %s',faultname));

% stairs(low,slip,PlotOpts_95);
% stairs(high,slip,PlotOpts_95);
% stairs(avg,slip,PlotOpts_avg);
% stairs(maxlikage,slip,PlotOpts_maxlik);
  

plot(high,slip,PlotOpts_95)
plot(low,slip,PlotOpts_95)
plot(avg,slip,PlotOpts_avg)
plot(maxlikage,slip,PlotOpts_maxlik)


