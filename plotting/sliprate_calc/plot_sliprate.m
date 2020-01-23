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
% You can also calculate and plot the slip rate (as a pdf) for a given time 
% in the model e.g. for the time-averaged slip rate for a given height
% (assinged by slip number). You can also plot the total slip rate
% (modelled scarp age/scarp height) and the elapsed time (time of the most
% recent EQ, avereaged between the models you have selected to show).
%
% The rates and timings will be plotted as a pdf (see ksdensity in matlab
% for more information on the nature of the distribution).
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
%write th number of models minus the burnin (e.g. 200k - 40k - write 160k)
%(This only applies to the sliphistory.out file, not TC's python script
%output - that is set in the setup.m file).
age_distribution = 160000;

%set file names

sliphistoryfile = '../../sliphistory.out';
%output from TC's python script
ratesfile = 'output_slip_averages_CAP.dat';

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

%slip rate and model numbers
 slips = load(ratesfile);
 timeS = slips(:,1);
 rateS = slips(:,2);
 modeldensity = slips(:,4) ./ slips(1,7);
 
% Maxlik calculation
% Sort Matrix based on likelihood and take the desired age distribution
[values, order] = sort(age(:,1));
age = age(order,:);
age = flipud(age);
age = age(1:age_distribution,4:end);

%Assign the scarp age and ET columns
scarpage = age(:,1);
ET = age(:,end-1);

%% calculate average slip rate for each model

 %calculate total slip, ignoring trench in mm
 slipT = sum(slip(1,1:end-1))*10;
 %calculate average slip rate for each slip history
 SR_avg = slipT ./ scarpage;

%Calculate mean values for slip rate, scarp age, and elapsed time
SR_mean = mean(SR_avg);
SA_mean = mean(scarpage);
ET_mean = mean(ET);
 
%calculate densities and quantiles - SRA
[DenseSRA,VariableSRA] = ksdensity(SR_avg);
QuantsSRA = quantile(SR_avg, [0.05 0.5 0.95]);

%calculate densities and quantiles - scarp age - SAA
[DenseSAA,VariableSAA] = ksdensity(scarpage);
QuantsSAA = quantile(scarpage, [0.05 0.5 0.95]);

%%calculate slip rate for a given height (e.g. the first 5 slips)
divider = (age(:,end-slipno)-age(:,end-slipno2))./2;

age_H = age(:,end-slipno)-divider(:,1);

slipH = sum(slip(1,end-slipno2:end-1))*10;

%slipH in metres
slipHm = slipH / 1000;

SR_H = slipH ./ age_H;
% 
% %remove weird really fast slip rates - if needed
% sortSRH = sort(SR_H);
% SR_H = sortSRH(1:end-200,:);

%calculate densities and quantiles - slip rate at height - SRH
[DenseSRH,VariableSRH] = ksdensity(SR_H);
QuantsSRH = quantile(SR_H, [0.05 0.5 0.95]);

%calculate densities and quantiles - slip age at height - SAH
[DenseSAH,VariableSAH] = ksdensity(age_H);
QuantsSAH = quantile(age_H, [0.05 0.5 0.95]);

%calculate densities and quantiles - elapsed time - ET
[DenseET,VariableET] = ksdensity(ET);
QuantsET = quantile(ET, [0.05 0.5 0.95]);

%maxlik calculations for plotting mean values
ET_max = age(2,end-1);
SA_max = age(2,1);
SR_max = slipT / SA_max;
%SR_H_max calculation
dividerM = (age(2,end-slipno)-age(2,end-slipno2))./2;
age_HM = age(2,end-slipno)-dividerM(:,1);
SR_H_max = slipH ./ age_HM;

%% plotting

%plot average slip rate distribution and SA
figure (54)

hold on

subplot(2,2,1);
area(VariableSRA,DenseSRA,'FaceColor',[0 0.4 0.4],'FaceAlpha',0.45)
hold on
y = get(gca,'YLim');
plot([QuantsSRA(1,2),QuantsSRA(1,2)],y,'b','LineWidth',1.2)
plot([SR_max,SR_max],y,'--k','LineWidth',1.2)
%plot([SRS,SRS],y,'Color',[0.91 0.41 0.17],'LineWidth',1.2)
xlabel('Average slip rate (mm/yr)')
%xlim([0 2])


subplot(2,2,2);
area(VariableSAA,DenseSAA,'FaceColor',[0.18 0.31 0.31],'FaceAlpha',0.45)
hold on
y = get(gca,'YLim');
plot([QuantsSAA(1,2),QuantsSAA(1,2)],y,'b','LineWidth',1.2)
plot([SA_max,SA_max],y,'--k','LineWidth',1.2)
%plot([SAS,SAS],y,'Color',[0.91 0.41 0.17],'LineWidth',1.2)
ax = gca;
ax.XRuler.Exponent = 0;
xlabel('Scarp age (years)')
%xlim([5000 24000])

%plot slip rate distribution for given height

subplot(2,2,3);
area(VariableSRH,DenseSRH,'FaceColor',[0 0.4 0.4],'FaceAlpha',0.45)
hold on
y = get(gca,'YLim');
plot([QuantsSRH(1,2),QuantsSRH(1,2)],y,'b','LineWidth',1.2)
plot([SR_H_max,SR_H_max],y,'--k','LineWidth',1.2)
%plot([SR_HS,SR_HS],y,'Color',[0.91 0.41 0.17],'LineWidth',1.2)
xlabel(sprintf('Slip rate, lower %.1f m (mm/yr)',slipHm))
%xlim([0 6])



subplot(2,2,4);
%area(VariableSAH,DenseSAH,'FaceColor',[0.18 0.31 0.31],'FaceAlpha',0.45)
area(VariableET,DenseET,'FaceColor',[0.18 0.31 0.31],'FaceAlpha',0.45)
hold on
y = get(gca,'YLim');
plot([QuantsET(1,2),QuantsET(1,2)],y,'b','LineWidth',1.2)
plot([ET_max,ET_max],y,'--k','LineWidth',1.2)
%plot([age_H_mean,age_H_mean],y,'b')
%plot([ETS,ETS],y,'Color',[0.91 0.41 0.17],'LineWidth',1.2)
xlabel('Elapsed time (years)')
%xlim([0 2500])


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




