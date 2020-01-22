%--------------------------------------------------------------------------
% Plotting function for output from MCMC_Flex_CPSS.m
% This produces a plot of the trace as well as the densities of each
% parameter as well as the statistical measurements
% 
% Run as: MCMCPlotTraces('../sliphistory.out',1) to include the entire 
% trace. If BurnIn = n > 1 the first samples n samples are ignore when 
% calculating the quantiles and the posterior distribution.
% The only inputs needed here are the prior used in the inference run.
% These will then be included as shaded area behind the pdfs
%--------------------------------------------------------------------------
function PlotMCMC = MMMCPlotTraces(TraceFileName,BurnIn)

    % Define the prior age settings used in the inference run. 
    % Mean and stddev of elapsed time (use -1 if uniform prior):
    MU_ET   = -1;
    SD_ET   = -1;
    % Mean and stddev of scarp age:
    MU_SA   = 15000;
    SD_SA   = 2500;

    %----------------------------------------------------------------------
    % NO FURTHER INPUT REQUIRED BEYOND THIS POINT
    %----------------------------------------------------------------------

    % Read the output from AgeOptim_MCMC_v3
    Trace = dlmread(TraceFileName);
    CPTrace = dlmread('../ncp.out');
    % Define plotlabels
    labvec = cell(5,1);
    labvec{1} = 'Elapsed Time';
    labvec{2} = 'Scarp Age';
    labvec{3} = 'Likelihood';
    labvec{4} = 'RMSW';
    labvec{5} = 'AICC';
    labvec{6} = '# Change points';
    
    % Activate and clear figure 1
    figure(1);
    clf;
    
    % Set up matrix to store all the quantiles
    QuantSummary = zeros(5,3);
    
    % Define vectors for prior parameters
    PriorMUVec = [MU_ET MU_SA -1 -1 -1]
    PriorSDVec = [SD_ET SD_SA -1 -1 -1]
    
    % Work out which parameters were estimated
    EstVec  = [5,6,7];
    for i = [2,1]
        if min(Trace(:,i)) < max(Trace(:,i))
            EstVec = [i,EstVec];
        end
    end
    nTraces = length(EstVec)  
    
    DataSize = size(Trace);
    ColumnNo = [(DataSize(2)-1) 4 1 2 3]
    
    % Run through all traces
    for p = 1:5
        % set subplot
        h(p) = subplot(6,3,3*(p-1)+1:3*(p-1)+2);
        % calculate densities and quantiles
        [Dens,Variable] = ksdensity(Trace(BurnIn+1:end,ColumnNo(p)));
        Quants = quantile(Trace(BurnIn+1:end,ColumnNo(p)), [0.05 0.5 0.95]);
        QuantSummary(p,1:3) = Quants;
        NoOfIts = size(Trace,1);
        
        % plot trace, 90% CI and median
        hold on;
        plot(Trace(2:end,ColumnNo(p)));set(gca,'FontSize',8);
        if p < 3
            plot([1 NoOfIts],[Quants(1) Quants(1)],'g')
            plot([1 NoOfIts],[Quants(2) Quants(2)],'r')
            plot([1 NoOfIts],[Quants(3) Quants(3)],'g')
        end
        axis([1 NoOfIts 0.99*min(Trace(2:end,ColumnNo(p))) 1.01*max(Trace(2:end,ColumnNo(p)))])
        
        if p==6
            plot(CPTrace)
        end
        hold off;
        xlabel('Iterations','FontSize',10);
        ylabel(labvec(p),'FontSize',10);
        
        % set subplot, and plot density, 90% CI and median
        h(10+p) = subplot(6,3,3*(p-1)+3);

        hold on;
        
        if p == 1
            area(Variable,normpdf(Variable,PriorMUVec(1),PriorSDVec(1)),...
                'FaceColor',[0.9, 0.9, 0.9],'EdgeColor','w')
        end
        if p == 2
            area(Variable,normpdf(Variable,PriorMUVec(2),PriorSDVec(2)),...
                'FaceColor',[0.9, 0.9, 0.9],'EdgeColor','w')
        end
        
        plot(Variable,Dens);set(gca,'FontSize',8);
        if p < 5
            plot([Quants(1) Quants(1)],[0 1.1*max(Dens)],'g')
            plot([Quants(2) Quants(2)],[0 1.1*max(Dens)],'r')
            plot([Quants(3) Quants(3)],[0 1.1*max(Dens)],'g')
        end
        

        axis([0.95*min(Trace(BurnIn+1:end,ColumnNo(p))) 1.05*max(Trace(BurnIn+1:end,ColumnNo(p))) 0 1.05*max(Dens)])
        ylabel('Density','FontSize',10)
        xlabel(labvec(p),'FontSize',10)
        
        
        
        hold off;
    end
    
    p = 6;
    % set subplot
    h(p) = subplot(6,3,3*(p-1)+1:3*(p-1)+2);
    % calculate densities and quantiles
    [Dens,Variable] = ksdensity(CPTrace(BurnIn:end));
    Quants = quantile(CPTrace(BurnIn:end), [0.05 0.5 0.95]);
    QuantSummary(p,1:3) = Quants;
    NoOfIts = size(Trace,1);
    
    % plot trace, 90% CI and median
    hold on;
    plot(CPTrace(1:end));set(gca,'FontSize',8);
    %         plot(Trace(:,p),'black');
    if p < 3
        plot([1 NoOfIts],[Quants(1) Quants(1)],'g')
        plot([1 NoOfIts],[Quants(2) Quants(2)],'r')
        plot([1 NoOfIts],[Quants(3) Quants(3)],'g')
    end
    min(CPTrace(1:end))
    max(CPTrace(1:end))
    axis([1 NoOfIts 0.99*min(CPTrace(1:end)) 1.01*max(CPTrace(1:end))])
    
    hold off;
    xlabel('Iterations','FontSize',10);
    ylabel(labvec(p),'FontSize',10);
    
    % set subplot, and plot density, 90% CI and median
    h(10+p) = subplot(6,3,3*(p-1)+3);
    hold on;
    plot(Variable,Dens);set(gca,'FontSize',8);
    %         plot(Variable,Dens,'black');
    if p < 5
        plot([Quants(1) Quants(1)],[0 1.1*max(Dens)],'g')
        plot([Quants(2) Quants(2)],[0 1.1*max(Dens)],'r')
        plot([Quants(3) Quants(3)],[0 1.1*max(Dens)],'g')
    end
    
    axis([0.95*min(CPTrace(BurnIn:end)) 1.05*max(CPTrace(BurnIn:end)) 0 1.05*max(Dens)])
    ylabel('Density','FontSize',10)
    xlabel(labvec(p),'FontSize',10)
    
    hold off;
    
    QuantSummary
    
    % count number of accepted changes
    Accepted = zeros(1,3);
    for i = 2:size(Trace,1)
        for j = 1:3
            if not(Trace(i,j) == Trace((i-1),j))
                Accepted(j) = Accepted(j) + 1;
            end
        end
    end
    
    NoAttempted = NoOfIts
    Accepted
    PropAccepted = sum(Accepted)/NoAttempted
           
% -------------------------------------------------------------------------
    % Generating plots of ages vs likelihood and RMSw
    figure(4)
    
    TraceSize = size(Trace)
    ColToPlot = [TraceSize(2)-1 4]
        
    for i = 1:2
        subplot(2,2,i)
        plot(Trace(:,ColToPlot(i)),Trace(:,1),'.')
        hold
        plot(Trace(end,ColToPlot(i)),Trace(end,1),'.','Color','r')
        hold off
        xlabel(labvec(i))
        ylabel('Likelihood')
    end
    for i = 1:2
        subplot(2,2,i+2)
        plot(Trace(:,ColToPlot(i)),Trace(:,2),'.')
        hold
        plot(Trace(end,ColToPlot(i)),Trace(end,2),'.','Color','r')
        hold off
        xlabel(labvec(i))
        ylabel('RMSw')
    end
    
end
