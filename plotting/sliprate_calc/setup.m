%%Setup for slip rate calculator and plotting
%
% Use this code to generate the input for TC's slip rate calculator in 
% python (slips.py). 
%
% The format for TC's code inserts 1 cm of slip at time = 0 and 1 cm of
% slip at time = scarp age. This allows for an effective rate to be
% calculated between the last earthquake and the modern day (and the
% oldest/second oldest event). This 'effective rate' would change if there
% was an earthquake tomorrow. You can change this by setting 'today' and
% 'startslip' to = 0.
%
% slips.py will produce an output file with the average slip rate binned at
% user-selected intervals, which can be plotted using plot_sliprate in this
% same folder.
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
%write th number of models minus the burnin (e.g. 200k - 40k - write 160k)
age_distribution = 10000;

%set file names

sliphistoryfile = '../../sliphistory.out';
%input to TC's python script
pyth_inp = 'slips_CAP.inp';

% set slip at 'today' and at the 'start' (e.g. oldest) to calculate an
% effective slip rate in TC's code for the start of your slip histories.

today = 1;

%% read in sliphistory and slip rate files and assign values

sliphistories = load(sliphistoryfile);
slip=sliphistories(1,4:end);
age=sliphistories(burnin:end,:); 
SizeOfFile = size(age);
maxmodel = SizeOfFile;

% Maxlik calculation
% Sort Matrix based on likelihood and take the desired age distribution
[values, order] = sort(age(:,1));
age = age(order,:);
age = flipud(age);
age = age(1:age_distribution,4:end);
 
% Format input file
 slips_format = [slip; age];
 slips_format(:,end) = today;
 dlmwrite(pyth_inp,slips_format,'delimiter',' ');





