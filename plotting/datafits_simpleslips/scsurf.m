function Ss = scsurf(Z,H,VARS)

% get required VARS:
dphi_dtheta = VARS.dphi_dtheta;
m_one_two_pi = VARS.m_one_two_pi;
rho_rock_dr_Lambda = VARS.scsurf.rho_rock_dr_Lambda;
sin_theta_m_cos_theta_theta_lt_B_theta_gt_C = ...
  VARS.sin_theta_m_cos_theta_theta_lt_B_theta_gt_C;
S_air = VARS.scsurf.S_air;

%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%
%--------------------------- scsurf.m -------------------------------------
%
% Calculates the scaling factor Ss for the exhumed samples
% as a function of:
%   Z = depth (cm) measured on the scarp (0 at surface, > 0 above),
%   H = height of the sarp (cm),
%   Lambda = the true attenuation length (g.cm-2) (for ex. 208 for neutrons), 
%   beta = scarp dip (degrees),
%   gamma = dip of upper eroded part of the scarp, above beta (degrees),
%   rho_rock = density (g.cm-3) of the rock.
%--------------------------------------------------------------------------

hz=H-Z; hz(hz==0) = 0.001; % to avoid NaNs

dr = ((exp(-(hz)*rho_rock_dr_Lambda)) .* sin_theta_m_cos_theta_theta_lt_B_theta_gt_C) * dphi_dtheta;

S_rock = (sum(dr(:)))*(m_one_two_pi) ;
Ss = S_air + S_rock ;

end
