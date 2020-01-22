%%Slip rate and slip rate variability calculator
%
% Use this code to plot the output from TC's slip rate calculator in 
% python (slips.py). 
%
% **This script needs some work!**
%
% The format for TC's code inserts 1 cm of slip at time = 0 and 1 cm of
% slip at time = scarp age. This allows for an effective rate to be
% calculated between the last earthquake and the modern day (and the
% oldest/second oldest event). This 'effective rate' would change if there
% was an earthquake tomorrow. You can change this by setting 'today' and
% 'startslip' to = 0.
%
% You can calculate the slip rate (as a pdf) for a given time in the
% model - e.g. for the average slip rate between 2 assigned slip numbers,
% for the average sliprate between the current day and a slip number, or
% for the average slip rate of the entired distribution (e.g. from today
% until the last slip). 
%
% The rate will be plotted as a pdf, which represents the distribution of
% rate calculated for the given slip number.
%
% 
%
% LG Jan 2020
clear
close all

%% Parameters 

%Fault name
faultname = 'Caporciano';

%set the burnin to remove the initial modelling phase before models settle
burnin = 40000;

%set the part of the distribution you would like to plot, e.g. for the top 
%1000 most likely slip histories, pick 1000. For the full distribution, 
%write 'end'
age_distribution = 10000;

%set file names

sliphistoryfile = '../../sliphistory.out';
%output from TC's python script
ratesfile = 'output_slip_averagesCAP.dat';

%set slip number for calculation. This counts back from the most recent
%event (which is slipno = 1), the penultimate event is slipno = 2, etc. The
%code will also output the average slip rate for the entire distribution,
%calculated by dividing the scarp age by the height of the scarp.

slipno = 6;
slipno2 = 5;

% set slip at 'today' and at the 'start' (e.g. oldest) to calculate an
% effective slip rate in TC's code for the start and end of your
% sliphistories.

today = 1;



%% read in sliphistory and slip rate files and assign values

sliphistory = load(sliphistoryfile);
slip=sliphistory(1,4:end);
age=sliphistory(burnin:end,4:end); 
scarpage = age(:,1);

%slip rate and model numbers
 slips = load(ratesfile);
 timeS = slips(:,1);
 rateS = slips(:,2);
 modeldensity = slips(:,4) ./ slips(1,7);

%% calculate average slip rate for each model

 %calculate total slip, ignoring trench in mm
 slipT = sum(slip(1,1:end-1))*10;
 %calculate average slip rate for each slip history
 SR_avg = slipT ./ scarpage;

 %calculate densities and quantiles - SRA
[DenseSRA,VariableSRA] = ksdensity(SR_avg);
QuantsSRA = quantile(SR_avg, [0.05 0.5 0.95]);

%calculate densities and quantiles - SAA
[DenseSAA,VariableSAA] = ksdensity(scarpage);
QuantsSAA = quantile(scarpage, [0.05 0.5 0.95]);

%%calculate slip rate for a given height (e.g. the first 5 slips)
divider = (age(:,end-slipno)-age(:,end-slipno2))./2;

age_H = age(:,end-slipno)-divider(:,1);

slipH = sum(slip(1,end-slipno2:end-1))*10;

SR_H = slipH ./ age_H;
% 
% %remove weird really fast slip rates - if needed
% sortSRH = sort(SR_H);
% SR_H = sortSRH(1:end-200,:);

%calculate densities and quantiles - SRH
[DenseSRH,VariableSRH] = ksdensity(SR_H);
QuantsSRH = quantile(SR_H, [0.05 0.5 0.95]);

%calculate densities and quantiles - SAH
[DenseSAH,VariableSAH] = ksdensity(age_H);
QuantsSAH = quantile(age_H, [0.05 0.5 0.95]);

%% plotting

%plot average slip rate distribution and SA
figure (54)

subplot(4,1,1);
area(VariableSRA,DenseSRA,'FaceColor',[0 0.4 0.4],'FaceAlpha',0.45)
hold on
xlim([0 10])
xlabel('Average slip rate (mm/yr)')

subplot(4,1,2);
area(VariableSAA,DenseSAA,'FaceColor',[0.18 0.31 0.31],'FaceAlpha',0.45)
xlim([0 24000])
xlabel('Scarp age, used in above (yrs)')

%plot slip rate distribution for given height

subplot(4,1,3);
area(VariableSRH,DenseSRH,'FaceColor',[0 0.4 0.4],'FaceAlpha',0.45)
xlim([0 10])
xlabel('Average slip rate, for lower 4.4 m of scarp (mm/yr)')

subplot(4,1,4);
area(VariableSAH,DenseSAH,'FaceColor',[0.18 0.31 0.31],'FaceAlpha',0.45)
xlim([0 3000])
xlabel('Time over which lower 4 m exhumed occurred (yrs)')

%plot slip rate over time
figure (55)

subplot(2,1,1);
area(timeS,rateS,'FaceColor',[0 0.4 0.4],'FaceAlpha',0.45)
hold on
%title('Average slip rate (mm/yr)')
xlabel('Time (yrs)')
ylabel('Slip rate (mm/yr)')

subplot(2,1,2);
area(timeS,modeldensity,'FaceColor',[0.18 0.31 0.31],'FaceAlpha',0.45)
ylabel('Density of models used in part 1')
xlabel('Time (yrs)')




