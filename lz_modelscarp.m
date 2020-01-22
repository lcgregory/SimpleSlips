%function Nf = modelscarp(data,coll,age,slip,preexp,EL,epsilon)
function [rmsw,rmswmean,lik,aicc,chi_square] = modelscarp(data,coll,age,slip,preexp,EL,epsilon,alpha,beta,gamma,rho_coll,Hfinal,VARS)

%-----------------------example--------------------------------------------
%     * first load the dataset :
% data = load('datarock.txt');
% coll = load('datacolluvium.txt');
% EL = load('datamagfield.txt');
%     * then, calculate :
% Nf = modelscarp(data,coll,age,slip,preexp,EL,epsilon);
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------
% modelscarp.m calculates the theoretical concentration in
% [36Cl]for samples caracterized by their chemical composition,
% position, thickness in the matrix 'data', after an earthquake sequence ;
% as well as RMS, AICC and I(%)
% plot of the [36Cl] dataset and of modeled [36Cl]
%
% data = [chemistry h thickness Cl_AMS sig_Cl_AMS] ;
% chemistry = data(:,1:n-4) ; chemistry : 62 colums
% h = data(:,n-3) ; samples positions (cm - !INTEGERS!) on scarp of dip beta
% if first h is at zero, then put h=1 to avoid NaNs.
% thickness = data(:,n-2) ; sample thickness when sampling (cm)
% Cl_AMS = data(:,n-1) ; [36Cl] sample concentration from AMS measurement
% sig_Cl_AMS = data(:,n) ; [36Cl] uncertainty on AMS measurement
%
% age = earthquakes (eq) ages, first is the oldest (yrs)
% slip = coseismic slip (cm) on fault scarp of dip beta
% preexp = pre-exposure duration before 1rst earthquake (yrs)
% (if samples from the buried part of the scarp were collected, put last
% eq at time zero, and put last 'slip' to the total height of collected samples)
%
% EL is a 4 columns matrix containing the epochs, associated time steps,
% and coefficients S_el,f (EL_f) and S_el,mu (EL_mu) which are function
% of elevation, latitude, longitude and intensity of Earth magnetic field
% at the study site (Dunai 2001; Pigati and Lifton, 2004; Lifton et al., 2005;
% Lifton et al., 2008 and references therein).
% For a constant Earth mag field use Stone 2000 and fix S_el,f and S_el,mu
% constants with time in the 'EL' file.
%       REMARK : Psi_Cl36_Ca_0 value from Stone et al., 1996 (48.8) was
% calculated for a constant mag field. This value must be changed below for
% each description of Earth mag field used (see clcoll.m and clrock.m)
%
% epsilon = erosion rate of scarp surface of dip beta (mm/yr)
%--------------------------------------------------------------------------
% 'f' stands for fast neutrons and 'mu' for slow muons
%--------------------------------------------------------------------------
% if pre-exp = 0, then comment the 'pre-exposition' part => part B
%--------------------------------------------------------------------------
%
% - A -
%
%------- ! PARAMETERS TO MODIFY FOR EACH SITE ! ---------------------------
%
% colluvial wedge dip alpha (degrees)
% alpha = 33 ; % MA3 : 30
% scarp dip beta (degrees)
% beta = 60 ; % MA3 : 45
% upper surface dip gamma (degrees)
% gamma = 33 ; % MA3 : 30
% present height of preserved scarp of dip beta at t = 0 (cm)
% Hfinal = 1100 ; % MA3 : 2000
% colluvial wedge mean density
% rho_coll = 1.5 ; % MA3 : 1.5
% rock (samples) mean density
rho_rock = VARS.rho_rock;
% Psi_Cl36_Ca_0 : Spallation production rate at surface of 40Ca
% ! depends on Earth Mag field description used !
%       Stone 2000: 48.8 +/- 3.5
%       Dunai 2001: 53.7 +/- 3.9
%       Pigati and Lifton 2004 (Desilets and Zreda, 2003): 53.1 +/- 3.8
%       Lifton et al., 2005: 59.4 +/- 4.3
%       Pigati and Lifton 2004 (Desilets et al., 2006): 54.7 +/- 4.0
%       Lifton et al., 2008: 58.9 +/- 4.3
Psi_Cl36_Ca_0 = 48.8 ;% (at of Cl36 /g of Ca per yr)
%
%--------------------------------------------------------------------------

%---------------------CONSTANTS--------------------------------------------
% Radioactive decay constant for 36Cl (a-1)
lambda36 = 2.303e-6 ;
%
% True attenuation length for fast neutron (g.cm-2)
Lambda = VARS.Lambda;
%--------------------------------------------------------------------------

%----------------EARTH MAG FIELD LOADING-----------------------------------
% Loading of Earth magnetic field variations from file 'EL'
if preexp == 1,
    EL(2,:)=EL(1,:); EL(2,1)=1; EL(2,2)=1;
end
if age(1) == 1,
    EL(2,:)=EL(1,:); EL(2,1)=1; EL(2,2)=1;
end
ti = EL(:,1) ; % time period (years)
it = EL(:,2) ; % time steps (years) - should be 100 yrs
EL_f = EL(:,3) ; % scaling factor for neutrons (S_el,f)
EL_mu = EL(:,4) ; % scaling factor for muons (S_el,mu)
%--------------------------------------------------------------------------

%---------------VARIABLES INITIALIZATION-----------------------------------
[m,n] = size(data) ; if n ~= 66, error('File data must have 66 columns'), end
nc = size(coll,2) ; if nc ~= 62, error('File coll must have 62 columns'), end
nel = size(EL,2) ; if nel ~= 4, error('File EL must have 4 columns'), end
%
N_eq = length(age) ; % number of earthquakes
%
R = sum(slip) ; % total cumulative slip
Rc = cumsum(slip) ;  Rc = [0 Rc] ;% slip added up after each earthquake
%
if Hfinal < R , error('Hfinal cannot be lower than cumulative slip R') , end
%
Hinit = Hfinal - R ; % initial height of the scarp during pre-exposure
%
%--------------------------------------------------------------------------

%--------------------SURFACE SCALING---------------------------------------
%               using scsurf.m for z>=0
% Calculates a scaling factor S_S(z>=0) every cm used for the samples at
% surface which is normalized by S_S(z=0) after in the calculation of production
% at surface (Parts B and C). This allows to take into account for the
% presence of upper part of dip gamma.
%
Zs = 0:1:R' ; % initialization of Zs ; one calculation point every cm
S_S = zeros(size(Zs)) ; % initialization of S_S (Surface Scaling)
%
for i = 1:length(Zs)    % loop on Zs
    a = scsurf(Zs(i),Hfinal,VARS) ;
    S_S(i) = a ;
end
%--------------------------------------------------------------------------

%----------------------DEPTH SCALING FOR NEUTRONS--------------------------
%       using scdepth.m, function of Hiseg (earthquakes) for z<=0
% Calculates a scaling factor S_D(z<=0) every 10 cm fitted by fitexp.m
% (S_D=so_f.exp(-z/Lambda_f). The derived so_f and Lambda_f depend
% on the height of the scarp of dip beta which grows after each earthquake
% (Hiseg = Hinitial + Rc(i) with earthquake i), so that so_f_d_iseg and
% Lambda_f_d_iseg are calculated iteratively.
% They are used later (parts B and C) to calculate the productions at depth
% which are then scaled to production at z=0 to derive a scaling factor
% function of z<=0.
%
so_f_diseg = zeros(N_eq,1) ; % initialization of so_f_d_iseg
Lambda_f_diseg = zeros(N_eq,1) ; % initialization of Lambda_f_d_iseg
%

Ziseg = 0:10:R' ; % one calculation point every 10 cm is sufficient
Ziseg = -Ziseg ; % negative because at depth

% calculate sc values:
Ziseg(Ziseg==0) = -0.001 ; % prevents NaNs
Zin_a = cell(size(Ziseg));
Zin_b = cell(size(Ziseg));
VARS.sc = cell(size(Ziseg));
for i = 1:length(Ziseg)
  Zin = Ziseg(i);
  Zin(Zin==0) = -0.001 ; % prevents NaNs
  Zin_a{i} = ((exp(Zin * VARS.scdepth.rho_coll_dv_Lambda)).*(VARS.sin_theta_m_cos_theta))*VARS.dphi_dtheta;
  Zin_a{i} = sum(Zin_a{i}(:));
  Zin_b{i} = ((exp(Zin * VARS.scdepth.rho_coll_da_Lambda)).*(VARS.sin_theta_m_cos_theta_theta_gt_B))*VARS.dphi_dtheta;
  Zin_b{i} = sum(Zin_b{i}(:));
  VARS.sc{i} = (Zin_a{i} + Zin_b{i})*(VARS.m_one_two_pi) ;
end

for is = 1:N_eq     % earthquake loop
    Hiseg = Hinit + Rc(is) ; % Height of exhumed scarp after each earthquake
    S_D_iseg = zeros(size(Ziseg)) ; % initialization of S_D_iseg
    for i = 1:length(Ziseg)     % loop on z
        VARS.scdepth.dr = exp(-(Hiseg - Ziseg(i))*VARS.scdepth.rho_rock_dr_Lambda) ; % H - Z > H
        S_D_iseg(i) = scdepth(i,VARS);
    end
    [dd,ee] = fitexp(-Ziseg*rho_coll,S_D_iseg,Lambda) ; % fit by fitexp.m
    so_f_diseg(is) = dd ; % constant so
    Lambda_f_diseg(is) = ee ; % attenuation length for neutron in direction z
end
%
% attenuation length perpendicular to colluvium surface after each
% earthquake (with H increasing after each earthquake):
Lambda_f_diseg = Lambda_f_diseg*sind(beta - alpha) ;

%---------------------------
% For beta infinite plane (used in B2 and C6):
Zbeta_inf = 0:10:1000'; Zbeta_inf = -Zbeta_inf ; % initialization
S_D_beta_inf = zeros(size(Zbeta_inf));
%

% calculate sc values:
Zbeta_inf(Zbeta_inf==0) = -0.001 ; % prevents NaNs
Zin_a = cell(size(Zbeta_inf));
Zin_b = cell(size(Zbeta_inf));
VARS.sc = cell(size(Zbeta_inf));
for i = 1:length(Zbeta_inf)
  Zin = Zbeta_inf(i);
  Zin_a{i} = ((exp(Zin * VARS.scdepth.rho_coll_dv_Lambda)).*(VARS.sin_theta_m_cos_theta))*VARS.dphi_dtheta;
  Zin_a{i} = sum(Zin_a{i}(:));
  Zin_b{i} = ((exp(Zin * VARS.scdepth.rho_coll_da_Lambda)).*(VARS.sin_theta_m_cos_theta_theta_gt_B))*VARS.dphi_dtheta;
  Zin_b{i} = sum(Zin_b{i}(:));
  VARS.sc{i} = (Zin_a{i} + Zin_b{i})*(VARS.m_one_two_pi) ;
end

for i = 1:length(Zbeta_inf)     % loop on z
        VARS.scdepth.dr = exp(-(2000 - Zbeta_inf(i))*VARS.scdepth.rho_rock_dr_Lambda) ; % H - Z > H
        S_D_beta_inf(i) = scdepth(i,VARS) ;
end
[so_f_beta_inf,Lambda_f_beta_inf] = fitexp(-Zbeta_inf*rho_coll,S_D_beta_inf,Lambda) ; % fit by fitexp.m
Lambda_f_beta_inf = Lambda_f_beta_inf*sind(beta - alpha) ; % attenuation perp. to colluvium surface
%--------------------------------------------------------------------------

%-------------------ROCK SCALING FOR NEUTRONS------------------------------
%        using scrock.m (attenuation in the direction of 'e')
%
e = 0:1:100' ; % e is in cm and perpendicular to scarp surface
Se = zeros(size(e)) ; % initialization of scaling Se
for i = 1:length(e)         % Loop on e
	Se(i) = scrock(e(i),VARS) ;
end
%
[so_f_e,Lambda_f_e] = fitexp(e*rho_rock,Se,Lambda) ; % exponential fit
%--------------------------------------------------------------------------

%---------------VARIABLES INITIALIZATION-----------------------------------
%
% h must be in cm and integers
h = data(:,n-3) ;  % initial positions of the samples at surface (cm)- integer
Z = (R - data(:,n-3))*rho_coll ; % initial depth of the samples (g.cm-2)

d = data ; % substitution of matrix data by matrix d
d(:,n-3) = Z ; % samples position along z
d(:,n-2) = data(:,n-2)*rho_rock ; % thickness converted in g.cm-2

slip_gcm2 = slip*rho_coll ; % coseismic slip in g.cm-2
sc = cumsum(slip_gcm2) ; % cumulative slip after each earthquake (g.cm-2)
sc0 = [0 sc] ;

% Positions along e initially (eo)
thick = data(:,n-2) ;
th2 = (thick/2)*rho_rock ; % 1/2 thickness converted in g.cm-2
eo = zeros(size(Z)) ;
for iseg = 1:N_eq
    eo(Z > sc0(iseg) & Z <= sc0(iseg + 1)) = epsilon*age(iseg)*0.1*rho_rock ; % in g.cm-2
end
eo(length(Z)) = epsilon*age(1)*0.1*rho_rock ;
eo = eo + th2 ; % we add the 1/2 thickness : sample position along e is given at the sample center
%--------------------------------------------------------------------------

%----- B ------------------------------------------------------------------
% comment pre-exposure part if pre-exp = 0, and uncomment the line below:
% N_in = zeros(size(Z)) ; Ni = zeros(size(Z)) ; Nf = zeros(size(Z)) ;
%-----------------------------PRE-EXPOSURE PROFILE-------------------------
%
% Calculation of [36Cl] concentration profile at the end of pre-exposure.
%
% initialization at 0
No = zeros(size(Z)) ; % No : initial concentration (here = zero) before pre-exposure
Ni = zeros(size(Z)) ; % Ni :
Nf = zeros(size(Z)) ; % Nf : final 36Cl concentration
%
tt = find(ti <= (age(1) + preexp) & ti > age(1)) ; % epoch index corresponding to pre-exposure
ip = it(tt) ; % corresponding intervals

% B1 - Loop - iteration on every samples
for j = 1:m

	dpj = d(j,:) ;
	dpj(n-3) = dpj(n-3)*sin((beta - alpha)*pi/180) ; % in the direction perpendicular to colluvium surface
	d0 = dpj ;
	d0(n-3) = 0 ;

% scaling is separated in two expressions: exp(-z/Lambda_f)*exp(-e/Lambda_f_e).

% in z direction : samples are at Z(j) = R - h(j)
% in e direction : samples are at position eo(j) = epsilon*T(j), where T(j) is the age of the earthquake which brings the sample j at surface.

	N_in = No(j) ; % initial concentration (here = zero)

        % calculate values for cl* functions ... dpj:
        VARS.chimie_dpj = dpj(1:(size(dpj, 2) - 4));
        VARS.chimie_dpj_44_53_chimie_dpj_58 = dpj(44) + dpj(45) + dpj(46) + dpj(47) + dpj(48) + dpj(49) + dpj(50) + dpj(51) + dpj(52) + dpj(53) + dpj(58);
        VARS.ppm44_dpj = dpj(44)*VARS.A_k(44)/(VARS.A_k(44) + 2*VARS.A_k(59));
        % d0:
        VARS.chimie_d0 = d0(1:(size(d0, 2) - 4));
        VARS.chimie_d0_44_53_chimie_d0_58 = d0(44) + d0(45) + d0(46) + d0(47) + d0(48) + d0(49) + d0(50) + d0(51) + d0(52) + d0(53) + d0(58);
        VARS.ppm44_d0 = d0(44)*VARS.A_k(44)/(VARS.A_k(44) + 2*VARS.A_k(59));

	% B2 - LOOP - iteration on time (ii) during pre-exposure
    for ii = 1:length(tt)

        % set to dpj values:
        VARS.chimie = VARS.chimie_dpj;
        VARS.chimie_44_53_chimie_58 = VARS.chimie_dpj_44_53_chimie_dpj_58;
        VARS.ppm44 = VARS.ppm44_dpj;

        [P_cosmo,P_rad] = clrock(d(j,:),eo(j),Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,VARS) ;
        % scaling at depth due to the presence of the colluvium: scoll=Pcoll(j)/Pcoll(z=0)
        P_coll = clcoll(coll,dpj,Lambda_f_diseg(1),so_f_diseg(1),EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,VARS) ;

        % switch to d0 values:
        VARS.chimie = VARS.chimie_d0;
        VARS.chimie_44_53_chimie_58 = VARS.chimie_d0_44_53_chimie_d0_58;
        VARS.ppm44 = VARS.ppm44_d0;

        P_zero = clcoll(coll,d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,VARS) ;
        scoll = P_coll/P_zero ;

        P_tot = P_rad + P_cosmo*scoll; % only P (Pcosmogenic) is scalled by scoll
        N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step
        N_in = N_out ;
    end

	Ni(j) = N_out ;

end
%--------------------------------------------------------------------------


%----- C ------------------------------------------------------------------
%-----------------------------SEISMIC PHASE--------------------------------
%
% -the term 'segment' is used for the samples exhumed by an earthquake-
% Calculation of [36Cl] profiles during seismic cycle.
% Separated in two stages :
%   * when samples are at depth and progressively rising because of earthquakes
%   (moving in the direction z with their position in direction e fixed)
%   * and when samples are brought to surface and only sustaining erosion
%   (moving along the direction e)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% FIRST EXHUMED SEGMENT is treated alone.
%
% variables initialization:
j1 = find(Z >= sc0(1) & Z <= sc0(2)) ; % samples from first exhumed segment
N1 = zeros(size(Z(j1))) ;
tt = find(ti <= age(1)) ; % epoch index more recent than first earthquake
ip = it(tt) ; % time intervals corresponding
%
% C1 - Loop - iteration on samples (k) from first exhumed segment
for k = 1:length(j1)

    djk = d(j1(k),:) ;
    hjk = h(j1(k)) ;    % position of sample k (cm)
    N_in = Ni(j1(k)) ;  % initial concentration is Ni, obtained after pre-exposure
    ejk = eo(j1(k)) ;   % initial position along e is eo(j1(k))

    % calculate values for clrock:
    VARS.chimie = djk(1:(size(djk, 2) - 4));
    VARS.chimie_44_53_chimie_58 = djk(44) + djk(45) + djk(46) + djk(47) + djk(48) + djk(49) + djk(50) + djk(51) + djk(52) + djk(53) + djk(58);
    VARS.ppm44 = djk(44)*VARS.A_k(44)/(VARS.A_k(44) + 2*VARS.A_k(59));

    % C2 - Loop - iteration on  time steps ii from t1 (= age eq1) to present
    for ii = 1:length(tt)
        [P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,VARS) ;

        scorr = S_S(1+hjk)/S_S(1) ;     % surface scaling factor (scorr)
        P_tot = P_rad + P_cosmo*scorr ;           % only Pcosmogenic is scaled with scorr
        N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step

        ejk = ejk - epsilon*ip(ii)*0.1*rho_rock ; % new position along e at each time step (g.cm-2)
        N_in = N_out ;
    end

    N1(k) = N_out ;

end

Nf(j1) = N1 ;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% ITERATION ON SEGMENTS 2 to N_eq
%
% C3 - Loop - iteration on each segment (from segment 2 to N_eq=number of eq)
for iseg = 2:N_eq

    j = find(Z > sc0(iseg) & Z <= sc0(iseg+1)) ; % index of samples from segment iseg
    z_j = Z(j) ; % initial depth along z of these samples (g.cm-2)
    N_new = zeros(size(z_j)) ;

    % C4 - Loop - iteration each sample from segment iseg
    for k = 1:length(j)

        ejk = eo(j(k)) ; % initial position along e is stil eo.
        djk = d(j(k),:) ;
        djk(n-3) = djk(n-3)*sind(beta - alpha) ;

        N_in = Ni(j(k)) ; %  initial concentration is Ni

        % C5 - Loop - iteration on previous earthquakes
        for l = 1:iseg-1
            ttt = find(ti <= age(l) & ti > age(l+1)) ; % epoch index
            ipp = it(ttt) ; % time intervals corresponding

            % depth (along z) are modified after each earthquake
            djk(n-3) = djk(n-3) - slip(l)*rho_coll*sind(beta - alpha) ;
			d0 = djk ;
			d0(n-3) = 0 ;

            % calculate values for cl* functions ... djp:
            VARS.chimie_djk = djk(1:(size(djk, 2) - 4));
            VARS.chimie_djk_44_53_chimie_djk_58 = djk(44) + djk(45) + djk(46) + djk(47) + djk(48) + djk(49) + djk(50) + djk(51) + djk(52) + djk(53) + djk(58);
            VARS.ppm44_djk = djk(44)*VARS.A_k(44)/(VARS.A_k(44) + 2*VARS.A_k(59));
            % d0:
            VARS.chimie_d0 = d0(1:(size(d0, 2) - 4));
            VARS.chimie_d0_44_53_chimie_d0_58 = d0(44) + d0(45) + d0(46) + d0(47) + d0(48) + d0(49) + d0(50) + d0(51) + d0(52) + d0(53) + d0(58);
            VARS.ppm44_d0 = d0(44)*VARS.A_k(44)/(VARS.A_k(44) + 2*VARS.A_k(59));

%------------------------------
            % C6 - DEPTH LOOP - iteration during BURIED PERIOD (T1 -> T(iseg-1))
%------------------------------
            for iii = 1:length(ttt)

                % set these to djk values:
                VARS.chimie = VARS.chimie_djk;
                VARS.chimie_44_53_chimie_58 = VARS.chimie_djk_44_53_chimie_djk_58;
                VARS.ppm44 = VARS.ppm44_djk;

            	[P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,VARS) ;
                % scaling at depth due to the presence of the colluvium: scoll=Pcoll(j)/Pcoll(z=0)
                P_coll = clcoll(coll,djk,Lambda_f_diseg(l+1),so_f_diseg(l+1),EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,VARS) ;

                % switch to d0 values:
                VARS.chimie = VARS.chimie_d0;
                VARS.chimie_44_53_chimie_58 = VARS.chimie_d0_44_53_chimie_d0_58;
                VARS.ppm44 = VARS.ppm44_d0;

                P_zero = clcoll(coll,d0,Lambda_f_beta_inf,so_f_beta_inf,EL_f(ttt(iii)),EL_mu(ttt(iii)),Psi_Cl36_Ca_0,rho_rock,rho_coll,VARS) ;
                scoll = P_coll/P_zero ;

                P_tot = P_rad + P_cosmo*scoll; % only P (Pcosmogenic) is scalled by scoll
                N_out = N_in + (P_tot - lambda36*N_in)*ipp(iii) ; % minus radioactive decrease during same time step
                N_in = N_out ;
            end

            N_in = N_out ;

        end

        N_in = N_out ;

        tt = find(ti <= age(iseg)) ; % epoch index more recent than earthquake iseg
        ip = it(tt) ; % time intervals corresponding
        djk = d(j(k),:) ;
        hjk = h(j(k)) ;

%------------------------------
            % C7 - SURFACE LOOP - iteration during EXHUMED PERIOD
%------------------------------

        % switch to djk values:
        VARS.chimie = djk(1:(size(djk, 2) - 4));
        VARS.chimie_44_53_chimie_58 = djk(44) + djk(45) + djk(46) + djk(47) + djk(48) + djk(49) + djk(50) + djk(51) + djk(52) + djk(53) + djk(58);
        VARS.ppm44 = djk(44)*VARS.A_k(44)/(VARS.A_k(44) + 2*VARS.A_k(59));

            for ii = 1:length(tt)
                [P_cosmo,P_rad] = clrock(djk,ejk,Lambda_f_e,so_f_e,EL_f(tt(ii)),EL_mu(tt(ii)),Psi_Cl36_Ca_0,rho_rock,VARS) ;

                scorr = S_S(1+hjk)/S_S(1) ; % surface scaling factor (scorr)
                P_tot = P_rad + P_cosmo*scorr ; % only Pcosmogenic is scaled with scorr
                N_out = N_in + (P_tot - lambda36*N_in)*ip(ii) ; % minus radioactive decrease during same time step

                ejk = ejk - epsilon*ip(ii)*0.1*rho_rock ; % new position along e at each time step (g.cm-2)
                N_in = N_out ;
            end

        N_new(k) = N_out ;

    end

    Nf(j) = N_new ;

end

cl36AMS = d(:,n-1) ; % sample concentration in [36Cl] measured by AMS
sig_cl36AMS = d(:,n) ; % uncertainty on [36Cl] AMS measurements

% RMSw (weighted least square) :
rmsw = ((cl36AMS - Nf)./sig_cl36AMS).^2 ;
rmswmean = mean(rmsw);
rmswmean = sqrt(rmswmean);
rmsw = sum(rmsw) ; % rmsw = sum(rmsw)/m ;
rmsw = sqrt(rmsw);

deviate = cl36AMS - Nf;
sqdeviate = deviate.^2;
lik     = exp(-(sum(sqdeviate))/(2*var(cl36AMS)));
% divisor = 1/((2*pi*(var(cl36AMS))^((size(cl36AMS,1))/2)))
% lik     = lik/((2*pi*(var(cl36AMS))^((size(cl36AMS,1))/2)))

% AICC (Akaike Information Criterion):
% if file data contains samples from the buried part of the scarp
% then, nb_param = 2*N_eq + 3 "-2" ;
% (and we add an earthquake at time = 0 and of slip = height of buried samples

nb_param = 2*N_eq + 3 ; % 2*N_eq + pre-exp + variance + erosion (epsilon)
hauteur = h;
if any(age == 0)
    nb_param = nb_param - 2 ;
    hauteur = h - slip(end) ;
end

aicc = ak(cl36AMS,Nf,nb_param);

% Chi_square
chi_square = ((cl36AMS - Nf)./sig_cl36AMS).^2 ;
chi_square = sum(chi_square) ;
chi_square = (1/(m - nb_param - 1))*chi_square;

% error bar coordinates
X = [cl36AMS(:)-sig_cl36AMS(:) cl36AMS(:)+sig_cl36AMS(:)] ;
X = X' ;
Y = [hauteur(:) hauteur(:)] ;
Y = Y' ;

% plot of calculated and measured values of 36Cl with associated error bars
plot(Nf,hauteur/100,'bo',cl36AMS,hauteur/100,'k.',X,Y/100,'k-') , hold on

hp = ones(N_eq,1)*get(gca,'XLim') ;
hp = hp' ;
vp = cumsum(fliplr(slip)) ;
vp = [vp' vp'] ;

if any(age == 0) % in case of samples coming from the buried part of the scarp
    vp = vp' - slip(end) ;
else
    vp = vp' ;
end

% plot of segments limits in z
plot(hp,vp/100,'b-') , hold off

% proportion of points where |calc - ams| <= sigma
i = find(abs(cl36AMS - Nf) <= sig_cl36AMS) ;
i = fix(100*length(i)/m) ;

% titles and legends
str_age = ['T (ka) = ',num2str(age/1000,'%5.1f')] ;
str_slip = ['; S (m) = ',num2str(slip/100,'%5.2f'),'; P (ka) = ',int2str(preexp/1000)] ;
str_geom = ['; a = ',int2str(alpha),'; b = ',int2str(beta),'; c = ',int2str(gamma)] ;
str_rho = ['; rrock = ',num2str(rho_rock),'; rcoll = ',num2str(rho_coll)] ;
str_rms = ['; RMSw = ',int2str(fix(rmsw)),'; Chi2 = ',int2str(chi_square),'; i = ',int2str(i),'; AICc = ',int2str(fix(aicc))] ;

%title([str_age str_slip str_geom str_rho str_rms])
xlabel('36Cl (at.g-1)')
% ylabel('Height on fault scarp (m)')
ylabel('Height (m)')
% axis([0.5e5 5e5 -4 12]) % axis for MA3 site (Magnola fault, Italy)

%----------- AICC FUNCTION ---------
function aicc = ak(measurements,calculations,K)

n = length(measurements) ;
aicc = sum((measurements - calculations).^2) ;
aicc = n*log(aicc/n) + (2*n*K)/(n - K - 1) ;
%-----------------------------------

%-lz--- AGE COMPARISON FUNCTION ----
