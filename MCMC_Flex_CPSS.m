%--------------------------------------------------------------------------
% This script is a front end to the modelscarp code. It uses a modified
% version of the original Schlagenhaug code (lz_modelscarp). The only
% difference being that a.o. the geometrical parameters are now defined in
% this script. There should be no need to change any parameters in
% lz_modelscarp.

% This script can only be run in Bayesian mode. In each iteration it either
% proposes, with some user-defined probability, a change to the elapsed
% time, the scarp age or the slip changepoints. Changes to the slip
% changepoints can be one of the following: insert a new slip changepoint,
% remove an existing slip changepoint or change the timing of an existing
% slip changepoint.

% Changes to the elapsed time and the total scarp age these are accepted or
% rejected using a standard Metropolis-Hastings algorithm. Changes to the
% changepoints are dealt with using  a reversible jump algorithm.

% Slip histories are generated using a constant of set and these are spaced
% in time in such a way that the resulting slip rate between change points
% is constant.

%--------------------------------------------------------------------------
function [RMSMat,AICCMat,ChiSqMat] = AgeOptim_MCMC

    %--------------------------------------------------------------------------
    % Preload data: make sure these files are in the correct location and
    % update for each fault
    data = load('datarockCAP.txt');
    coll = load('datacolluviumCAP.txt');
    EL = load('datamagfieldCAP_SfTV.txt');

    % Preset some parameters
    epsilon = 0;
    preexp = 100;

    % colluvial wedge dip alpha (degrees)
    alpha = 22.0 ;
    % scarp dip beta (degrees)
    beta = 48.5 ;
    % upper surface dip gamma (degrees)
    gamma = 28.7 ;

    % colluvial wedge mean density
    rho_coll = 1.5;

    %--------------------------------------------------------------------------
    % Present height of preserved scarp of dip beta at t = 0 (cm), does not
    % include trench depth
    ScarpHeight = 820;
    % Depth of trench below scarp (cm)
    TrenchDepth = 45;

    % Present height of preserved scarp of dip beta at t = 0 (cm)
    Hfinal = ScarpHeight + TrenchDepth;

    % A sequence of slip events is calculated given the scarp height and
    % changepoint heights above. Offsets (cm) are drawn from a uniform
    % distribution (Ave-error,Ave+error) until the entire scarp is covered.
    % Set error to 0 for fixed offsets
    AveOffset       = 102.5;
    OffsetError     = 0; % Values other than zero are not compatible with
                         % this version of the code

    % Error on inter event times, this assumes inter event times vary
    % randomly by a factor of about +- TimeError (0.0 <= TimeError < 1).
    % Note that when TimeError = 0 the inter event times are fixed between
    % change points. Do not use TimeError >= 1 as this causes inter to be
    % zero or negative
    TimeError = 0.0;     % Values other than zero are not compatible with
                         % this version of the code

    % Estimate elapsed time (2 = normal prior, 1 = uniform, 0 = don't
    % estimate)
    EstET = 1;

    % Estimate scarp age (2 = normal prior, 1 = uniform. No option for 0)
    % LG 10/18
    EstSA = 2;

    %----------------------------------------------------------------------
    % Parameters for bayesian inference
    % Start new run or append to existing (1 = new, 0 = append)
    StartNew = 1;

    % File to which to write output, note that if restarting run this
    % filename has to match a pre-existing output file
    FileToUse   = 'sliphistory.out';
    FileToRead  = 'sliphistory.out';
    CPFileToUse = 'ncp.out';
    TCPFileToUse = 'tcp.out';%HG 2017 creates time of change point file
    MCPFileToUse = 'mcp.out';%HG 2017 creates matrix of change point file

    % Number of iterations to run estimation
    nIts = 200000;

    % Initial ages
    InitET          = 500
    InitAge         = 15000

    % Standard deviations for proposal distributions
    PropSigAge      = 400;

    % Prior for Elapsed time, if ETSD == 0.0 elapsed time is not estimated
    ETMean          = 1000
    ETSD            = 1000

    % Prior for scarp age. If ETSA = 1, then this is not used
    ScarpAgeMean    = 15000
    ScarpAgeSD      = 2500

    % Prior for scarp age uniform distribution. If ETSA = 1, this is set
    % such that all SA within the minimum and maximum have equal
    % probability, and everything SA outside of the min and max have 0
    % probability. USE WITH CAUTION, set wide bounds, and this will
    % probably never converge %LG 10/18
    ScarpAgeMin     = 8000;
    ScarpAgeMax     = 24000;

    % Define probabilities for changing ET, SA or all other ages
    % respectively, it is generally not necessary to change these.....
    if EstET == 0
        ChangeProbs = [0 1/10 1];
    else
        ChangeProbs = [1/10 2/10 1];
    end

    % --- set up various variables here:

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

    %----------------------------------------------------------------------
    % NO FURTHER INPUT REQUIRED BEYOND THIS POINT
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    % Initialising the MCMC run

    % If restarting a run read parameters included in estimates from
    % preexisting output. Note that this doesn't check that the other
    % parameters are identical, so some care has to be taken here.
    if StartNew == 0
        PreData = dlmread(FileToRead,' ');
        Pre_nCP = dlmread(CPFileToUse,' ');
        size(PreData);

        CurrLik     = PreData(end,1);
        CurrRMSW    = PreData(end,2);
        CurrAICC    = PreData(end,3);
        CurrnCP     = Pre_nCP(end);

        slip        = PreData(1,4:end);
        SlipSize    = size(slip);
        nOffsets    = SlipSize(2)-1;
        age         = PreData(end,4:end)
        CurrAge     = age;
        PropAge     = CurrAge;

        LikTrace    = PreData(:,1);
        RMSWTrace   = PreData(:,2);
        AICCTrace   = PreData(:,3);
        nAccepted   = 0;
        for i = 2:size(LikTrace,1)
            if not(LikTrace(i) == LikTrace(i-1))
                nAccepted = nAccepted + 1;
            end
        end
    else
        % Calculate slip vector
        [slip]      = LinearCreateSlip(TrenchDepth,ScarpHeight,AveOffset,OffsetError,Hfinal);
        SlipSize    = size(slip);
        nOffsets = SlipSize(2)-1;
        % Calculate age vector
        age = CreateAges(InitAge,InitET,TimeError,nOffsets);
        CurrAge     = age;
        PropAge     = CurrAge;
        CurrLik     = 1e-50;
        CurrRMSW    = NaN;
        CurrAICC    = NaN;
        LikTrace    = NaN;
        RMSWTrace   = NaN;
        AICCTrace   = NaN;
        nAccepted   = 0;
        CurrnCP     = 0;
        PropnCP     = CurrnCP;
    end

    if EstET == 2
        PriorET     = normpdf(age(end-1),ETMean,ETSD);
    else
        PriorET     = 1;
    end

    if EstSA == 2               %LG 10/18
        PriorScarpAge   = normpdf(age(1),ScarpAgeMean,ScarpAgeSD);
    else
        PriorScarpAge   = pdf('Uniform',age(1),ScarpAgeMin,ScarpAgeMax);        %LG 10/18
    end

    SlipSize    = size(slip);
    nOffsets    = SlipSize(2)-1;

    % Create a vectors that stores an index for each slip event except the
    % first and the last two. Also create a vector that stores the 1 and
    % one but last offset. These vectors will be modified in the RJMC step:
    % When a change point is added, one of the offsets is moved from the
    % Offsets vector to the ChangePoints vector. The reverse happens when
    % a changepoint is removed. Changes in timing do not affect these
    % vectors
    if (StartNew == 0)
        [Offsets,ChangePoints,nChPoints] = FindOffsetsFromOld(CurrAge);
        PropOffsets  = Offsets;
        PropCPs      = ChangePoints;
        PropnCP      = nChPoints;
        CurrnCP      = PropnCP;
    else
        Offsets      = 2:SlipSize(2)-2;
        PropOffsets  = Offsets;
        ChangePoints = [1 SlipSize(2)-1];
        PropCPs      = ChangePoints;
        nChPoints    = 0;
    end

    % Calculate initial model and associated measures of fit
    if nIts == 0.0
        figure(10);
    end

    [rmsw,rmswmean,lik,aicc,chi_square] = lz_modelscarp(data,coll,CurrAge,slip,preexp,EL,epsilon,alpha,beta,gamma,rho_coll,Hfinal,VARS);
    CurrLik     = lik
    CurrRMSW    = rmsw
    CurrAICC    = aicc;

    % Write first two lines of output file
    if StartNew == 1
        dlmwrite(FileToUse,[NaN, NaN, NaN, slip],' ');
        dlmwrite(FileToUse,[lik, rmsw, aicc, CurrAge],'-append','delimiter',' ');
        dlmwrite('ncp.out',CurrnCP);
        dlmwrite('tcp.out',PropCPs); %HG 2017
        %dlmwrite('mcp.out',
    end

    % Run MCMC parameter estimation

    % Set random number generator
    rng(1,'twister');

    % Initialise counters to keep track of type of changes made
    CountET = 0;
    CountSA = 0;
    CountAdd = 0;
    CountRemove = 0;
    CountShift = 0;

    % Initialize PropProbRatio
    PropProbRatio = 1;

    % Run through all iterations
    for it=1:nIts
        it

        % Pick a random number between 0 and 1
        RandNo = rand;

        if RandNo < ChangeProbs(1)
            % Make a change to elapsed time and increment counter
            %LG changed - did say 'scarp age' but it is actually refering
            %to ET
            PropAge(end-1) = normrnd(PropAge(end-1), PropSigAge);
            CountET = CountET + 1;
            % Set proposal probability ratio
            PropProbRatio   = 1;
        elseif RandNo < ChangeProbs(2)
            % or make a change to the scarp age and increment counter
            %LG changed - did say 'elapsed time' but it is actually refering
            %to scarp age
            PropAge(1) = normrnd(PropAge(1), PropSigAge);
            CountSA = CountSA + 1;
            % Set proposal probability ratio
            PropProbRatio   = 1;
        else
            % or make a change to the changepoints and increment counter
            % TEMP define ModProbs
            ModProbs = [1/10 2/10 1];
            % Pick a random number between 0 and 1
            RandNo = rand;

            if RandNo < ModProbs(1)
                % Add a change point
                if numel(Offsets) > 1
                   % Choose index of 1 entry from Offsets
                    OffToChange = Offsets(randi(numel(Offsets)));
                    PropCPs     = sort([ChangePoints,OffToChange]);
                    PropOffsets = Offsets(Offsets~=OffToChange);
                    PropnCP   = PropnCP + 1;
                    % Increment counter
                    CountAdd = CountAdd + 1;
                    % Find index position of new entry in changepoint vector
                    CPIndex         = find(PropCPs==OffToChange);
                    % Find indices of new age and age before and after in age
                    % vector
                    AgeIndex        = PropCPs(CPIndex);
                    AgeIndexBefore  = PropCPs(CPIndex-1);
                    AgeIndexAfter   = PropCPs(CPIndex+1);
                    % Find corresponding times
                    CPAgeBefore = PropAge(AgeIndexBefore);
                    CPAgeAfter  = PropAge(AgeIndexAfter);
                    % Calculate time interval
                    TimeInt  = CPAgeBefore - CPAgeAfter;

                    % Work out time for new changepoint
                    PropAge(AgeIndex) = PropAge(AgeIndexBefore) - rand*TimeInt;

                    % Set proposal probability ratio
                    PropProbRatio = (ModProbs(2)-ModProbs(1)) / ModProbs(1);
                    PropProbRatio = PropProbRatio * TimeInt / (PropAge(1) - PropAge(end-1));
                    ProbProbRatio = PropProbRatio * (numel(Offsets)+1)/nChPoints;

                end

                % Pick random time between them
            elseif RandNo < ModProbs(2)
                % Remove a changepoint
                if numel(ChangePoints) > 2
                    % Choose index of 1 entry from Offsets
                    CPToRemove = ChangePoints(randi(numel(ChangePoints)-2)+1);

                    % Calculate time interval, needed for proposal
                    % probability ratio
                    % Find index position of new entry in changepoint vector
                    CPIndex         = find(PropCPs==CPToRemove);
                    % Find indices of new age and age before and after in age
                    % vector
                    AgeIndex        = PropCPs(CPIndex);
                    AgeIndexBefore  = PropCPs(CPIndex-1);
                    AgeIndexAfter   = PropCPs(CPIndex+1);
                    % Find corresponding times
                    CPAgeBefore = PropAge(AgeIndexBefore);
                    CPAgeAfter  = PropAge(AgeIndexAfter);
                    % Calculate time interval
                    TimeInt  = CPAgeBefore - CPAgeAfter;

                    PropOffsets = sort([Offsets,CPToRemove]);
                    PropCPs = PropCPs(ChangePoints~=CPToRemove);
                    PropnCP = PropnCP - 1;
                    % Increment counter
                    CountRemove = CountRemove + 1;

                    % Times will be updated later when recalculating age vector

                    % Set proposal probability ratio
                    PropProbRatio = ModProbs(1) / (ModProbs(2)-ModProbs(1));
                    PropProbRatio = PropProbRatio * (PropAge(1) - PropAge(end-1)) / TimeInt;
                    ProbProbRatio = PropProbRatio * (nChPoints+1)/numel(Offsets);
%                     ProbProbRatio = 1


                end

            else
                % Change timing
                if numel(ChangePoints) > 2
                    % Choose index of 1 entry from Offsets
                    CPToChange = PropCPs(randi(numel(PropCPs)-2)+1);
                    % Increment counter
                    CountShift = CountShift +1;
                    % Find index position of entry in changepoint vector
                    CPIndex         = find(PropCPs==CPToChange);
                    % Find indices of new age and age before and after in age
                    % vector
                    AgeIndex        = PropCPs(CPIndex);
                    AgeIndexBefore  = PropCPs(CPIndex-1);
                    AgeIndexAfter   = PropCPs(CPIndex+1);
                    % Find corresponding times
                    CPAgeBefore = PropAge(AgeIndexBefore);
                    CPAgeAfter  = PropAge(AgeIndexAfter);
                    % Calculate time interval
                    TimeInt  = CPAgeBefore - CPAgeAfter;
                    % Work out time for new changepoint
                    PropAge(AgeIndex) = PropAge(AgeIndexBefore) - rand*TimeInt;

                    % Set proposal probability ratio
                    PropProbRatio = 1;
                end

            end
        end
        % Recalculate ages
        PropAge = CreateNewAges(PropAge,PropCPs);

        % Report the types of changes made
        [CountET  CountSA  CountAdd CountRemove CountShift];

        % Need to check that ages are positive and in correct order
        AgeDiffs = PropAge(1:nOffsets)-PropAge(2:nOffsets+1);

        if min(AgeDiffs) >= 0.0 %LG 2018 - corrected bug where AgeDiffs was PropAge in orginal
            age = PropAge;

            % run schlagenhauf code
            [rmsw,rmswmean,lik,aicc,chi_square] = lz_modelscarp(data,coll,PropAge,slip,preexp,EL,epsilon,alpha,beta,gamma,rho_coll,Hfinal,VARS);

            % Calculate prior for ET
            if EstET == 2
                PriorPropET     = normpdf(PropAge(end-1),ETMean,ETSD);
            else
                PriorPropET     = 1;
            end

            % Calculate prior for scarp age %LG 10/18
            if EstSA == 2
                PriorPropScarpAge   = normpdf(PropAge(1),ScarpAgeMean,ScarpAgeSD);
            else
                PriorPropScarpAge   = pdf('Uniform',PropAge(1),ScarpAgeMin,ScarpAgeMax); %LG 10/18
            end

            % Calculate acceptance ratio in stages
            AcceptRatio = lik/CurrLik;
            AcceptRatio = AcceptRatio*PriorPropET/PriorET;
            AcceptRatio = AcceptRatio*PriorPropScarpAge/PriorScarpAge;
            AcceptRatio = AcceptRatio*PropProbRatio;
            AcceptProb  = min(1.0,AcceptRatio);
        else

            AcceptProb = 0.0;
        end

        % Work out whether to accept or reject proposed change
        RandNo = rand;
        if RandNo < AcceptProb
            CurrAge         = PropAge;
            Offsets         = PropOffsets;
            ChangePoints    = PropCPs;
            CurrnCP         = PropnCP;
            CurrLik         = lik;
            PriorET         = PriorPropET;
            PriorScarpAge   = PriorPropScarpAge;
            CurrRMSW        = rmsw;
            CurrAICC        = aicc;
            nAccepted       = nAccepted + 1;
        else
            PropAge         = CurrAge;
            PropOffsets     = Offsets;
            PropCPs         = ChangePoints;
            PropnCP         = CurrnCP;
        end

        % Work out proportion accepted
        PropAccepted = (nAccepted-1)/it;

        tcp = CurrAge(ChangePoints); % HG 2017

        % Write output files
        dlmwrite(FileToUse,[CurrLik, CurrRMSW, CurrAICC, CurrAge],'-append','delimiter',' ');
        dlmwrite('ncp.out',CurrnCP,'-append');
        dlmwrite('tcp_index.out',ChangePoints,'-append');
        dlmwrite('tcp.out',tcp,'-append'); %HG 2017 - Writes time of change point to tcp.out
        dlmwrite('mcp.out',ChangePoints); % HG 2017 - Writes matrix of change point to mcp.out
        dlmwrite('PropAccepted.out',PropAccepted,'-append'); %LG - writes proportion of models accepted

    end

    figure(10);
    PlotCumSlip(age,slip);

end

%--------------------------------------------------------------------------
% This function works out the sequence of slip events. This is currently
% only done once but could be run in every iteration if offset is not
% fixed. This requires change to the MCMC routine.
function [slip] = LinearCreateSlip(TrenchDepth,ScarpHeight,AveOffset,OffsetError,Hfinal)

    % Calculation of slip vector working from the bottom.
    CumHeight   = TrenchDepth;
    slip        = [TrenchDepth];

    % Add events after second change height
    while CumHeight < Hfinal
        Offset      = AveOffset + OffsetError * (2*rand - 1);
        slip        = [Offset,slip];
        CumHeight   = CumHeight + Offset;
%         nEvents3    = nEvents3 + 1;
    end

    % Add one additional event to third part of history
    Offset      = AveOffset + OffsetError * (2*rand - 1);
    slip        = [Offset,slip];
    CumHeight   = CumHeight + Offset;

    % Need to adjust first element as total slip shouldn't exceed the scarp
    % height
    slip(1)     = Hfinal - sum(slip(2:end))- 0.001;

    if slip(1) <= 0
        slip = slip(2:end);
        slip(1)     = Hfinal - sum(slip(2:end))- 0.001;
    end
end

%--------------------------------------------------------------------------
% This function works out the age vector by generating a random number for
% each event in the slip vector these are scaled so that the change point
% ages and scarp age fit exactly.
function [age] = CreateAges(MaxAge,ElapsedTime,TimeError,nOffsets)

    % Calculate age vector

    % Define age vector, first age = 0, second age = elapsed time, rest
    % of ages equally spaced on total age range. Note that ages are
    % stored in reverse order.
    NoOfEvents  = nOffsets+1;
    AgeInc      = (MaxAge-ElapsedTime)/(NoOfEvents-2);

    % Generate inter event times with random error, note that these
    % are scaled so that they will add up to exactly the change age
    RVec       = rand((NoOfEvents-2),1);
    ScaledRVec = 1- TimeError + TimeError*RVec/mean(RVec);
    TIncVec    = AgeInc * ScaledRVec;

    age = [ElapsedTime,0];
    for i = 1:length(TIncVec)
        age = [(age(1)+TIncVec(i)),age];
    end
end

%--------------------------------------------------------------------------
% Find offset, changepoints and number of changepoints from given
% age-vector. This is slightly complicated by the fact that the slip
% histories are not stored with very high precision and thus the
% inter-event times are slightly affected by rounding errors
function [Offsets,ChangePoints,nChPoints] = FindOffsetsFromOld(CurrAge);

    % find number of slips
    nSlips = length(CurrAge)-1;

    % initialise vectors storing indices of change points and offsets
    ChangePoints = [1 nSlips];
    Offsets = 2:(nSlips-1);
    nChPoints = 0;
    % calculate inter-event times
    Diffs = CurrAge(1:end-3)-CurrAge(2:end-2);

    % run through all events
    for j = 1:nSlips-3
        if abs(Diffs(j)-Diffs(j+1)) > 1
            Offsets         = Offsets(Offsets~=(j+1));
            ChangePoints    = sort([ChangePoints, (j+1)]);
            nChPoints       = nChPoints + 1;
        end
    end
end

%--------------------------------------------------------------------------
% This function works out the age vector by from the vector of changepoiint
% indices. Sections between change points are assumed to have constant slip
% rate and thus ages are spaced out equally between them
function [age] = CreateNewAges(age,ChangePoints)

    % Find how many changes there are and run through all sections
    DimCP   = size(ChangePoints);
    for i = 2:DimCP(2)
        % Find index of changepoint after and calculate no of times to
        % change and new time increment
        CPBefore    = ChangePoints(i-1);
        CPAfter     = ChangePoints(i);
        nChanges    = CPAfter - CPBefore;
        TInc = (age(CPAfter) - age(CPBefore))/(nChanges);
        % Recalculate all ages
        for ti = CPBefore:CPAfter-2
            age(ti+1) = age(ti) + TInc;
        end
    end
end

%--------------------------------------------------------------------------
% This function generates a plot of the slip history. Only called from the
% end of the MCMC run.
function PlotCum = PlotCumSlip(age,slip)
    % Plot slip history
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
    age = [age(1),age];
    CumSlip	= [(CumSlip(1)+slip(1)),CumSlip];

    figure(11);
    stairs(age,CumSlip);
    xlabel('Time (yrs)');
    ylabel('slip (cm)');
end







