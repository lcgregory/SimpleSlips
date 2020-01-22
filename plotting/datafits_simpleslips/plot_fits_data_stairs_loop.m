%% script to plot fits to the data in order of liklihood and the corresponding 
%  slip history stairs
%
% Requires sliphistory.out file (in the above directory) 
% User must input the site specific details for the lz_modelscarp code here
%
% This file should be run in a folder with datafits_modelscarp
%
% LCG Oct 2018 + modification from RR on the MCMC/SimpleSlips code
%
% Tip: Type 'opengl hardware' into matlab before making a plot to have
% better resolution figures
% 
% Once you have run this script once for a site, you can just run the later
% plotting sections with the same ages.out and fits.out (or comment out the
% loop)

clear
close all
%% Set up for plotting fit to data

faultname = 'Caporciano';

% Input file names (from main SimpleSlips file) and set site parameters

data_file = '../../datarockCAP.txt';
coll_file = '../../datacolluviumCAP.txt';
EL_file = '../../datamagfieldCAP_SfTV.txt';

slip_file = '../../sliphistory.out';

% colluvial wedge dip alpha (degrees)
alpha = 22.0 ;
% scarp dip beta (degrees)
beta = 48.5 ;
% upper surface dip gamma (degrees)
gamma = 28.7 ;

% colluvial wedge mean density
rho_coll = 1.5;

% Present height of preserved scarp of dip beta at t = 0 (cm), does not
% include trench depth
ScarpHeight = 820;
% Depth of trench below scarp (cm)
TrenchDepth = 45;

% Present height of preserved scarp of dip beta at t = 0 (cm)
Hfinal = ScarpHeight + TrenchDepth;

% Other parameters in the Schlagenhauf code
epsilon = 0;
preexp = 100;

% Set up the datafits parameters
% choose the steps to take in running the iteration
steps = 1000;

% Set burnin
burnin = 40000;

% delete existing file - comment out if you are appending (e.g. have to
% stop a run)
 delete 'ages.out';
 delete 'fits.out';

%% No further input required below this point 

% Output file names
FitsFile = 'fits.out';
AgesFile = 'ages.out';

% Load sliphistory file
SlipHistories = dlmread(slip_file);

%Schlagenhauf requirements/lz_modelscarp requirements
data = load(data_file);
coll = load(coll_file);
EL = load(EL_file);


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
    
%% Set up variables here (from Richard Rigby's modification):

 % init VARS struct:
    VARS = [];

    % -- lz_modelscarp vars:
    Lambda = 208;
    rho_rock = 2.7;
    VARS.Lambda = Lambda;
    VARS.rho_rock = rho_rock;

    % -- cl*.m values:

    % A_k = atomic mass of element k (g.mol-1)
    A_k = [74.9 137.33 9.01218 209.0 112.4 140.1 58.9332 51.996 132.9054 63.5] ;
    A_k = [A_k 162.5 167.3 152.0 69.7 157.25 72.6 178.5 164.9 114.8 138.9] ;
    A_k = [A_k 175.0 95.94 92.9 144.2 58.7 207.2 140.9 85.4678 121.8 150.4] ;
    A_k = [A_k 118.7 87.62 180.9 158.9 232 168.9 238.029 50.9 183.8 88.9] ;
    A_k = [A_k 173.0 65.4 91.22 28.0855 26.98154 55.847 54.938 24.305 40.08 22.98977] ;
    A_k = [A_k 39.0983 47.9 30.97376 10.81 6.941 1.0079 32.06 12.011 15.9994 15.9994 35.453] ;

    % Num_k = Atomic number of element k
    Num_k = [33 56 4 83 48 58 27 24 55 29] ;
    Num_k = [Num_k 66 68 63 31 64 32 72 67 49 57] ;
    Num_k = [Num_k 71 42 41 60 28 82 59 37 51 62] ;
    Num_k = [Num_k 50 38 73 65 90 69 92 23 74 39] ;
    Num_k = [Num_k 70 30 40 14 13 26 25 12 20 11] ;
    Num_k = [Num_k 19 22 15 5 3 1 16 6 8 8 17] ;

    % Xi_k = average log-decrement of energy loss per collision for element k
    Xi_k = [0 0 0 0 0 0 0 0.038 0 0] ;
    Xi_k = [Xi_k 0 0 0 0 0.013 0 0 0 0 0] ;
    Xi_k = [Xi_k 0 0 0 0 0 0 0 0 0 0.013] ;
    Xi_k = [Xi_k 0 0 0 0 0 0 0 0 0 0] ;
    Xi_k = [Xi_k 0 0 0 0.07 0.072 0.035 0.036 0.08 0.049 0.084] ;
    Xi_k = [Xi_k 0.05 0.041 0 0.174 0.264 1 0 0.158 0.12 0.12 0.055] ;

    % sigma_sc_k = neutron scattering x-section of element k (barns)

    sigma_sc_k = [0 0 0 0 0 0 0 3.38 0 0] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 0 172 0 0 0 0 0] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 0 0 0 0 0 0 38] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 0 0 0 0 0 0 0] ;
    sigma_sc_k = [sigma_sc_k 0 0 0 2.04 1.41 11.35 2.2 3.42 2.93 3.025] ;
    sigma_sc_k = [sigma_sc_k 2.04 4.09 5 4.27 0.95 20.5 0 4.74 3.76 3.76 15.8] ;

    % sigma_th_k = thermal neutron absorbtion x-section of element k (barns)
    sigma_th_k = [0 0 0 0 0 0 0 3.1 0 0] ;
    sigma_th_k = [sigma_th_k 0 0 0 0 41560 0 0 0 0 0] ;
    sigma_th_k = [sigma_th_k 0 0 0 0 0 0 0 0 0 9640] ;
    sigma_th_k = [sigma_th_k 0 0 0 0 0 0 0 0 0 0] ;
    sigma_th_k = [sigma_th_k 0 0 0 0.17 0.23 2.56 13.3 0.063 0.43 0.53] ;
    sigma_th_k = [sigma_th_k 2.15 6.1 0.2 767 70.5 0.33 0 0.0034 0.0002 0.0002 33.5] ;

    % I_a_k = dilute resonance integral for absorption of epithermal neutrons by element k (barns)
    I_a_k = [0 0 0 0 0 0 0 1.6 0 0] ;
    I_a_k = [I_a_k 0 0 0 0 390 0 0 0 0 0] ;
    I_a_k = [I_a_k 0 0 0 0 0 0 0 0 0 1400] ;
    I_a_k = [I_a_k 0 0 0 0 0 0 0 0 0 0] ;
    I_a_k = [I_a_k 0 0 0 0.127 0.17 1.39 14 0.038 0.235 0.311] ;
    I_a_k = [I_a_k 1 3.1 0 1722 0 0 0 0.0016 0.0004 0.0004 13.7] ;

    % f_d_k = proportion of muons stopped in element k that are captured by the nucleus
    f_d_k = [0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0 0 0 0 0 0 0] ;
    f_d_k = [f_d_k 0 0 0 0.671 0.582 0.906 0 0.538 0.864 0.432] ;
    f_d_k = [f_d_k 0.83 0 0 0 0 0 0 0.09 0.223 0.223 0] ;

    % Y_n = average neutron yield per captured muon
    Y_n = [0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_n = [Y_n 0 0 0 0.86 1.26 1.125 0 0.6 0.75 1] ;
    Y_n = [Y_n 1.25 0 0 0 0 0 0 0.76 0.8 0.8 0] ;

    % S_i = mass stopping power (MeV/(g.cm-2))
    S_i = [0 0 0.000529 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0 0 0 0 0 0 0] ;
    S_i = [S_i 0 0 0 0.000454 0.000444 0.000351 0 0.000461 0.000428 0.000456] ;
    S_i = [S_i 0.000414 0.000375 0.000433 0.000527 0.000548 0 0.000439 0.000561 0.000527 0.000527 0] ;

    % Y_U_n = neutron yield (n/an/g/ppm de U)
    Y_U_n = [0 0 265 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_U_n = [Y_U_n 0 0 0 0.69 5.1 0.19 0 5.8 0 14.5] ;
    Y_U_n = [Y_U_n 0.45 0 0 62.3 21.1 0 0 0.45 0.23 0.23 0] ;

    % Y_TH_n = neutron yield (n/an/g/ppm de Th)
    Y_Th_n = [0 0 91.2 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0 0 0 0 0 0 0] ;
    Y_Th_n = [Y_Th_n 0 0 0 0.335 2.6 0.205 0 2.6 0 6.8] ;
    Y_Th_n = [Y_Th_n 0.305 0 0 19.2 9.6 0 0 0.18 0.079 0.079 0] ;

    VARS.A_k = A_k;
    VARS.Num_k = Num_k;
    VARS.Xi_k = Xi_k;
    VARS.sigma_sc_k = sigma_sc_k;
    VARS.sigma_th_k = sigma_th_k;
    VARS.I_a_k = I_a_k;
    VARS.f_d_k = f_d_k;
    VARS.Y_n = Y_n;
    VARS.S_i = S_i;
    VARS.Y_U_n = Y_U_n;
    VARS.Y_Th_n = Y_Th_n;

    VARS.coll_44_53_coll_58 = sum([coll(44:53) coll(58)]);

    % -- sc*.m vars:

    [theta,phi] = meshgrid(0:90,0:180) ;
    theta = theta*pi/180 ;
    phi = phi*pi/180 ;
    dphi = pi/180 ;
    dtheta = pi/180 ;

    % Lal exponent:
    m = 2.3;
    m_one = 2.3 + 1;

    two_pi = pi * 2;
    pi_180 = pi / 180;
    m_one_two_pi = m_one / two_pi;

    VARS.m_one_two_pi = m_one_two_pi;

    cos_theta = cos(theta);
    cos_theta_sin_phi = cos(theta).*sin(phi);
    sin_theta = sin(theta);
    sin_theta_m = sin_theta.^m;
    sin_theta_m_cos_theta = (sin_theta_m).*cos_theta;
    phi_pi = phi + pi;
    sin_phi_pi = sin(phi_pi);
    cos_theta_sin_phi_pi = cos_theta.*sin(phi_pi);
    sin_phi = sin(phi);
    dphi_dtheta = dphi * dtheta;
    sin_theta_m_cos_theta_dphi_theta = sin_theta_m_cos_theta*dphi_dtheta;
    sum_sin_theta_m_cos_theta_dphi_theta = sum(sin_theta_m_cos_theta_dphi_theta(:));

    VARS.sin_theta_m_cos_theta = sin_theta_m_cos_theta;
    VARS.dphi_dtheta = dphi_dtheta;

    alpha_pi_180 = alpha * pi_180;
    beta_pi_180 = beta * pi_180;
    gamma_pi_180 = gamma * pi_180;
    cos_beta = cos(beta_pi_180);
    sin_beta = sin(beta_pi_180);
    tan_beta = tan(beta_pi_180);
    tan_gamma = tan(gamma_pi_180);
    cos_gamma = cos(gamma_pi_180);
    sin_gamma = sin(gamma_pi_180);
    cos_alpha = cos(alpha_pi_180);
    sin_alpha = sin(alpha_pi_180);
    sin_beta_alpha = sin(beta_pi_180 - alpha_pi_180);
    sin_beta_gamma = sin(beta_pi_180 - gamma_pi_180);

    B = atan(tan_beta.*sin_phi);
    C = atan(tan_gamma.*sin_phi);

    theta_lt_B = theta < B;
    theta_gt_B = theta > B;
    theta_gt_C = theta > C;
    sin_theta_m_cos_theta_theta_gt_B = ...
      sin_theta_m_cos_theta.*theta_gt_B;
    sin_theta_m_cos_theta_theta_lt_B_theta_gt_C = ...
      sin_theta_m_cos_theta.*theta_lt_B.*theta_gt_C;

    VARS.sin_theta_m_cos_theta_theta_gt_B = sin_theta_m_cos_theta_theta_gt_B;
    VARS.sin_theta_m_cos_theta_theta_lt_B_theta_gt_C = sin_theta_m_cos_theta_theta_lt_B_theta_gt_C;

    scdepth = [];
    scdepth.dv = abs((sin_beta_alpha) ./ ...
                          (sin_theta.*cos_alpha - ...
                           sin_alpha.*cos_theta.*sin_phi_pi));
    scdepth.da = abs((sin_beta_alpha) ./ ...
                          (sin_theta.*cos_alpha - ...
                           sin_alpha.*cos_theta_sin_phi));
    scdepth.rho_coll_dv_Lambda = rho_coll * scdepth.dv / Lambda;
    scdepth.rho_coll_da_Lambda = rho_coll * scdepth.da / Lambda;
    scdepth.rho_rock_dr_Lambda = rho_rock * scdepth.da / Lambda;

    VARS.scdepth = scdepth;

    scsurf = [];
    scsurf.dv = sum_sin_theta_m_cos_theta_dphi_theta;
    scsurf.da = sin_theta_m_cos_theta_theta_gt_B*dphi_dtheta;
    scsurf.da = sum(scsurf.da(:));
    scsurf.dr = abs((sin_beta_gamma) ./ ...
                          (sin_theta.*cos_gamma - ...
                           sin_gamma.*cos_theta_sin_phi));
    scsurf.rho_rock_dr_Lambda = rho_rock * scsurf.dr / Lambda;
    scsurf.S_air = (scsurf.da + scsurf.dv)*(m_one_two_pi) ;

    VARS.scsurf = scsurf;

    scrock = [];
    scrock.da = abs(1 ./ (sin_theta.*cos_beta - ...
                               sin_beta.*cos_theta_sin_phi));
    scrock.dv = abs(1 ./ (sin_theta.*cos_beta - ...
                               sin_beta.*cos_theta_sin_phi_pi));
    scrock.rho_rock_da_Lambda = rho_rock * scrock.da / Lambda;
    scrock.rho_rock_dv_Lambda = rho_rock * scrock.dv / Lambda;

    VARS.scrock = scrock;

    % ---
    
%% Start of Loop

for array_number = 0:steps:maxmodel;
array_number
    
% Set age array from sorted slip histories

maxlik = SizeOfFile(1) - array_number;
age = SortedSlipHistories(maxlik,4:end);

% Run Schlagenhauf code for age and slip array

datafits_modelscarp(data,coll,age,slip,preexp,EL,epsilon,alpha,beta,gamma,rho_coll,Hfinal,VARS);

%Save fits and slips
 fits = ans';
 dlmwrite(FitsFile,[array_number,fits],'-append','delimiter',' ');
 dlmwrite(AgesFile,[array_number,age],'-append','delimiter',' ');

end

%% Plot 36Cl data and fits
cl36 = data(:,65);
cl36sigma = data(:,66);
height = data(:,63);
height = height - TrenchDepth;
unc = [cl36-cl36sigma,cl36+cl36sigma];

cl36Fits = load(FitsFile);
NumberOfModels = size(cl36Fits);


% Plot figure in order such that most likely fits turn up on top.
figure(28)
hold on
cc = parula(NumberOfModels(1));
cc = flipud(cc);
cl36Fits = flipud(cl36Fits);
for k = 1:NumberOfModels(1)
plot(cl36Fits(k,2:end),height,'color',cc(k,:));
end

title(sprintf('Measured and modelled ^{36}Cl: %s', faultname));
xlabel('^{36}Cl concentration (at/g)');
ylabel('Sample height (cm)');

% Circles filled in
plot(cl36,height,'o','color','k','MarkerFaceColor','w','DisplayName','Data')
% Open circles
% plot(cl36,height,'o','color','b','DisplayName','Data')
% legend
unc = unc';
height = [height,height];
height = height';
plot(unc,height,'color','k')


%% Plot slip histories that were used in the comparison

agesSorted = load(AgesFile);
ages = load(AgesFile);
NumberOfAges = size(ages);

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
    
    % Figure
    figure(32);
    hold on;
    
    ax = gca;
    ax.XRuler.Exponent = 0;
    
    xlabel('Time (yrs)');
    ylabel('slip (cm)');
    title(sprintf('Modelled slip histories for ^{36}Cl comparison: %s',faultname));
    
    % Plot slip histories in color scale according to order in file
    % (likelihood) - flipped so the 'best fits' are on top
    
     cc = parula(NumberOfAges(1)); 
     cc = flipud(cc);
     agesSorted = flipud(agesSorted);
     
     for i = 1:NumberOfAges(1)
         ages = agesSorted(i,2:end);
         ages = [ages(1),ages];
         stairs(ages,CumSlip,'color',cc(i,:));
     end
     
 
%% Plot most likely slip history on top - not currently coded in but could 
%  be added

% age_plot = [age(1),age];
% stairs(age_plot,CumSlip,'red','LineWidth',1.8);
% 
% ylim([0 maxheight]) %new for plotting all at same size
% xlim([0 23000])

% Add another slip history plotted on top of histoplot manually, using same
% slips as in maxlik file: 

% ageC = [15500 13562 11625 9687 7750 5812 3875 1937 0];
% ageC = [ageC(1),ageC];
% stairs(ageC,CumSlipM,'blue','LineWidth',1.8);
