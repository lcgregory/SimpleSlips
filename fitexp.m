function [so_f,Lambda_f] = fitexp(z,s,Lambda)

% [so_f,Lambda_f] = fitexp(z,s,Lambda)
%
%%--------------------------------------------------------------------------
% Schlagenhauf A., Gaudemer Y., Benedetti L., Manighetti I., Palumbo L.,
% Schimmelpfennig I., Finkel R., Pou K.
% G.J.Int., 2010
%-------------------------- ? ---------------------------------------------

% Calculates the coefficient Lambda_f
% of the exponential fit of the scaling s(z) calculated by 
% scdepth.m, scsurf.m and scrock.m :
% s(z = 0).exp(-z/Lambda)
%
% Lambda_f has the SAME units as z and Lambda (m, cm or g.cm-2)
%

% valeur initiale du paramètre Lambda
start_point = rand(1) ;

% le modèle est décrit par la fonction expfun
model = @expfun ;

options = optimset('MaxIter',1e+6) ;

% recherche des meilleurs paramètres (moindres carrés)
estimates = fminsearch(model,start_point,options) ;

so_f = s(1) ;
k = estimates(1) ;
Lambda_f = Lambda/k ;

% détail de la fonction expfun
    function [sse,expz] = expfun(params)
        % so = params(1) ;
        k = params(1) ;
        expz = s(1)*exp(-z*k/Lambda) ;
        misfit = expz - s ;
        sse = sum(misfit.^2) ;
    end

end