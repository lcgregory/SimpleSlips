readme file for plotting results from SimpleSlips

This folder contains several plotting scripts that have been developed by Leo Zijerveld (Bergen), Laura Gregory (Leeds), and Huw Goodall (Leeds). Laura has developed many of them, and she is not very good at Matlab. Apologies for the sometimes clunky coding.

You must first run MCMC_Flex_CPSS.m for your site. We recommend 100-500k iterations (~200k is standard for what we've done so far). See documentation in the code readme and Cowie et al., 2017, DOI: 10.1038/srep44858 

Each plotting script or function can run as a standalone routine. You need to first open the script and edit it for your site and what you want to plot. There are further instructions within each script. 

Take care with what you are plotting - double check that it is what you expected and don't over-interpret the results (e.g. your most likely slip history is not very representative of the potential exposure histories your site has experienced, but it is informative to see a pattern of slip that fits the data well). We have tried to develop plots that reflect the simplicity of our modelling.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_stairs.m
This shows the most likely slip histories plotted as individual lines, rather than as a density of lines (as in histo_plot). You can choose how many models to plot in the parula colour scale, and how many models to plot in grey behind it. For example, you may want to show the most likely 1000 slip histories in colour, on top of the 10,000 most likely histories in grey. It takes a lot of time to plot this many lines, so more than 10,000 will run slowly, but it is possible.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SimpleSlips_PlotTraces.m
This is a function that shows the trace of the model, the density of each parameter, and statistical measurements for each parameter, calculated excluding the burnin that you specifiy in running the function. Convergence is indicated by parameters returning to the same range in values on the lefthand plots.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Confidence interval
Both of these scripts plot the confidence interval of your distribution (minus a burnin). These lines could be plotted on top of other plots (e.g. the stairs.m) if you want to show both some of the detail and the range of the full accepted distribution. 

%smooth_confidence.m
This plots bounds on the distribution of the model as a high and low percentile. It is binned such that the percentile is calculated for age at each step in the slip. This code plots a smooth distribution in age vs time. The mean is also plotted (default red) and the most likely slip history (blue). 

%stepped_confidence.m
The same as the smooth confidence script, but the distribution is plotted as steps rather than smooth lines. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datafits_simpleslips

% plot_fits_data_stairs_single.m
Here you can plot the fit to the data of one model (e.g. your most likely model) and the slip vs. time for that exposure history. This is a useful plot for showing how the 'best' model does or does not fit the data.

% plot_fits_data_stairs_loop.m
Here you can loop through your full model distribution. Because the MCMC code does not output the fits of the data for each model, this steps through the attempted models to show a selection of the fit to the data (and the models that have been plotted). This gives you a sense of the range of fits that have been accepted, and the variability of the distribution. You should always check how your models fit the data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
histo_plot
%plot_slip_2D_hist.m 
This plots a 2D histogram of the models accepted by MCMC_Flex_CPSS and should be run in the same directory as hist2d.m. The histogram reflects the slip v time parameter space covered by the modelling, but does not show the detailed slip history. There are tradeoffs between scarp age and slip rate changes that may not be reflected in this plot. Plot stairs shows the tradeoff in the latter relationship better.

The second plot that comes up is slip rate vs time. This is more controversial (LG - might not include it in published version) because it is complicated. Each step back in time is the *average* slip rate, NOT the slip rate for that time interval. E.g. the point at 14 kyr is the slip rate averaged from 0-14kyr (for CAP ~0.5 mm/yr) not the slip rate AT that time in the modelling (closer to 0). This plot is also a 2D histogram, and thus reflects the distribution of the time-averaged slip rate of the accepted models.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sliprate_calc

*To be used with caution. Please read Gregory et al. (submitted), Goodall et al. (submitted), and Gregory et al. (in prep) to understand what assumptions are important for this plot. *

% slips.py
This is a python script written by Tim Craig (University of Leeds) to calculate the mean slip rate in bins in time (decided by the user), which can then be plotted as slip rate vs. time (see plot_sliprate.m). The mean is taken from a suite of models input by the user - see figure in the supplement of Gregory et al. (submitted) for a graphical representation. You have to run it in python, which is installed on most linux machines. It requires a formatted version of the sliphistory.out file, which is produced from:

% setup.m
Setup for slips.py

%plot_sliprate.m
This script takes the output from slips.py and plots sliprate as a function of time or slips in the model. Because the MCMC code proposes slip history models by only changing time, not slip magnitude, you can simply pick which portion of the slip rate you would like to plot as a pdf. You can also plot sliprate vs time, or the average slip rate over the length of the model.

The output from slips.py can be plotted in other methods (e.g. GMT was used to plot figure 3 in Gregory et al., submitted).



% % % % % % % % % % % % % % % % % % % % %
events.m

This hasn't been written yet. I plan to write a script that plots histograms of the 'events' on a time vs. density plot

