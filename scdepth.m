function sd = scdepth(i,VARS)

% get required VARS:
dr = VARS.scdepth.dr;
dphi_dtheta = VARS.dphi_dtheta;
m_one_two_pi = VARS.m_one_two_pi;
sin_theta_m_cos_theta_theta_lt_B_theta_gt_C = ...
  VARS.sin_theta_m_cos_theta_theta_lt_B_theta_gt_C;
sc = VARS.sc{i};

% sd = scdepth(Z,H,Lambda,alpha,beta,gamma,rho_rock,rho_coll)
%
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%
%--------------------------- scdepth.m ------------------------------------
%
% Calculates the scaling factor sd for the buried samples
% as a function of:
%   Z = depth (cm) measured on the scarp (0 at surface, < 0 underneath),
%   H = height of the scarp (cm),
%   Lambda = the true attenuation length (g.cm-2) (for ex. 208 for neutrons),
%   alpha = colluvium dip (degrees) ;
%   beta = scarp dip (degrees),
%   gamma = dip of upper eroded part of the scarp, above beta (degrees),
%   rho_rock = density (g.cm-3) of the rock,
%   rho_coll = density (g.cm-3) of the colluvium.
%--------------------------------------------------------------------------

dr = (dr .* sin_theta_m_cos_theta_theta_lt_B_theta_gt_C) * dphi_dtheta;
sr = (sum(dr(:))) * m_one_two_pi;
sd = sc + sr ;

end
