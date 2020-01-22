function sr = scrock(h,VARS)

% get required VARS:
da = VARS.scrock.da;
dphi_dtheta = VARS.dphi_dtheta;
m_one_two_pi = VARS.m_one_two_pi;
rho_rock_da_Lambda = VARS.scrock.rho_rock_da_Lambda;
rho_rock_dv_Lambda = VARS.scrock.rho_rock_dv_Lambda;
sin_theta_m_cos_theta = VARS.sin_theta_m_cos_theta;
sin_theta_m_cos_theta_theta_gt_B = ...
  VARS.sin_theta_m_cos_theta_theta_gt_B;

% sr = scrock(h,lambda,beta,rho_rock)
%
%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%
%----------------------------- scrock.m -----------------------------------
%
% attenuation in the direction of e (perpendicular to fault scarp of dip beta)
% as a function of:
%   h = position in direction e of the sample
%   Lambda = true attenuation length (g.cm-2) (for ex. 208 for neutrons),
%   beta = scarp dip (degrees)
%   rho_rock = density (g.cm-3) of the rock.
%--------------------------------------------------------------------------

h(h==0) = 0.001 ; % prevents NaNs 

da = ((exp(-h*rho_rock_da_Lambda)) .* sin_theta_m_cos_theta_theta_gt_B) * dphi_dtheta;
sa = (sum(da(:))) * m_one_two_pi;

dv = ((exp(-h*rho_rock_dv_Lambda)) .* sin_theta_m_cos_theta) * dphi_dtheta;
sv = (sum(dv(:))) * m_one_two_pi;

sr = sa + sv ;

end
