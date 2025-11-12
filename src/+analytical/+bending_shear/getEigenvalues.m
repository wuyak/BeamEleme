% function cheb_rep = getTimoshenkoEigenvalues(sys, BC)
% function for finding the eigenvalues
% it returns a chebfun represenation that can be qeried for roots.
% The algorithm uses exact, numerically stable expressions as described
%  in the manuscript:
% Exact and Numerically Stable Expressions for Euler-Bernoulli and
% Timoshenko Beam Modes, by Firas A. Khasawneh and Daniel J. Segalman
% The user also has the option of using the ill-conditioned conventional
% expressions.
% If you use this package, please cite the above paper, as well as the code
% repo:
%
% Inputs:
% sys: a structure which defines the system. It includes:
%   beam parameters:
%   sys.par.A: cross-sectional area of the beam
%   sys.par.I: mass moment of area of the beam
%   sys.par.L: length of the beam
%   %
%   non-dimensional parameters:
%   sys.par.alpha: the value of non-dimensional parameter alpha=A*L^2/I
%   sys.par.beta: the value of non-dimensional parameter beta=((A*G*ksq)/(Y*I))*L^2
%   sys.par.gamma: the value of the non-dimensional parameter
%                  gamma=beta/alpha=G*ksq/Y
%   sys.par.T: non-dimensional time scaling factor T_T= L*sqrt(rho/(G*ksq))
%   %
%   algorithm options:
%   sys.numStable: logic flag, if true, use numerically stable expressions,
%                  if false, use conventional (ill-conditioned) expressions
%   sys.MaxMode: the suprememum of the natural frequency to consider
%   sys.tol: the tolerance to use when comparing expressions to zero
%   optional parameters:
%   sys.lambda_scale: (optional) a parameter used for scaling the eigvals
%   sys.norm: (optional) how to normalize the mode shapes. The options are
%           'mass' for mass orthonormalization, or 'unscaled' for no scaling.
%   %
% BC: the boundary conditions of the beam. Options are not case sensitive,
% and the include:
%   BC = "PP": pinned-pinned (checked)
%   BC = "CF": clamped-free (checked)
%   BC = "CP": clamped-pinned (checked)
%   BC = "PF": pinned-free (checked)
%   BC = "CC": clamped-clamped (checked)
%   BC = "FF": free-free (checked)
%   BC = "PR": pinned-roller (checked)
%   BC = "CR": clamped-roller (checked)
%   BC = "RR": roller-roller (checked)
%   BC = "RF": roller-free (checked)
%  %
% Output:
% chep_rep: a 3x1 cell array of chebfunction representations of the
% characterisitc function. The cells contain the chebfun represenation for
% lambda < alpha, lambda=alpha, and lambda > alpha, respectively. The
% resulting represenatations can be quereis for roots using the roots()
% command.
%
function cheb_rep = getTimoshenkoEigenvalues(sys, BC)

% Add chebfun to path if needed
if ~exist('chebfun', 'file')
    chebfun_path = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'chebfun');
    if exist(chebfun_path, 'dir')
        addpath(genpath(chebfun_path));
    else
        error('chebfun not found. Please install from https://www.chebfun.org/ or check path: %s', chebfun_path);
    end
end

% get the scaling for lambda, if any
if ~isfield(sys, 'lambda_scale')
    x_scale = 1;
else
    x_scale = sys.lambda_scale;
end

% x*x_scale is lambda
BC = upper(BC);
tol = sys.tol;  % tolerance to check equality to zero

% extract needed parameters
alphaN = sys.par.alpha;
gammaN = sys.par.gamma;
% betaN = sys.par.beta;
maxMode = sys.maxMode;
numStable = sys.numStable;
wSqrt = sqrt(alphaN*(1+gammaN));

% get the needed useful params
deltaN = @(x) 1 - 4*gammaN/(1+gammaN).^2 .* (1-alphaN./x);
muN = @(x) sqrt(0.5*x*(1+gammaN)*(sqrt(deltaN(x))-1));
wN = @(x) sqrt(0.5*x*(1+gammaN)*(sqrt(deltaN(x))+1));
thetaN = @(x) sqrt(0.5*x*(1+gammaN)*(1-sqrt(deltaN(x))));
psiN = @(x) atan2(1-exp(-2*muN(x)), 1+exp(-2*muN(x)));
gN = @(x) sqrt(2)*exp(-muN(x))./sqrt(1+exp(-4*muN(x)));

range_LT = [tol, sys.par.alpha-tol]/x_scale;
range_GT = [sys.par.alpha+tol, max(maxMode, 1.2*sys.par.alpha)]/x_scale;

switch BC
    case "CC"    % fixed-fixed
        % obtain the corresponding eigenvalues
        if ~numStable  % this is the conventional expression
            eigVal_LT  = chebfun (@(x) sinh(muN(x*x_scale)).*sin(wN(x*x_scale)).*(((wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2))./(muN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)))-((muN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2))./(wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2))))+2*cosh(muN(x*x_scale)).*cos(wN(x*x_scale))-2, range_LT, 'splitting', 'on');
        else  % this is the stabilized expression
            eigVal_LT = chebfun(@(x) sin(psiN(x*x_scale)).*sin(wN(x*x_scale)).*((wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2))./(muN(x*x_scale)*(x*x_scale-wN(x*x_scale).^2))-((muN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2))./(wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2))))+ 2*cos(psiN(x*x_scale)).*cos(wN(x*x_scale))-2*gN(x*x_scale), ...
                range_LT, 'splitting', 'on');
        end
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(alphaN * sin(wSqrt) + 2 * gammaN*alphaN/wSqrt * (1-cos(wSqrt))) < tol;
        eigVal_EQ = eq_check * alphaN;
        eigVal_GT = chebfun (@(x) 2*(x*x_scale - thetaN(x*x_scale).^2) .* (x*x_scale-wN(x*x_scale).^2).*(cos(wN(x*x_scale)).*cos(thetaN(x*x_scale))-1)./(wN(x*x_scale).*thetaN(x*x_scale)) + sin(wN(x*x_scale)).*sin(thetaN(x*x_scale)).*((x*x_scale)^2.*wN(x*x_scale).^2 + (x*x_scale).^2*thetaN(x*x_scale).^2 - 4 * x*x_scale .* wN(x*x_scale).^2.*thetaN(x*x_scale).^2 + wN(x*x_scale).^4.*thetaN(x*x_scale).^2+wN(x*x_scale).^2.*thetaN(x*x_scale).^4)./(wN(x*x_scale).* thetaN(x*x_scale)).^2, ...
            range_GT, 'splitting', 'on');
        %----------------------------
    case "CF"    % fixed-free (cantilever beam)
        % obtain the corresponding eigenvalues
        if ~numStable  % this is the conventional expression
            eigVal_LT  = chebfun(@(x) ((x*x_scale + muN(x*x_scale).^2)./(x*x_scale - wN(x*x_scale).^2)+(x*x_scale - wN(x*x_scale).^2)./(x*x_scale + muN(x*x_scale).^2)).*cosh(muN(x*x_scale)).*cos(wN(x*x_scale)) + (wN(x*x_scale)./muN(x*x_scale) - muN(x*x_scale)./wN(x*x_scale)).*sinh(muN(x*x_scale)).*sin(wN(x*x_scale))-2, range_LT, 'splitting', 'on');
        else  % this is the stabilized expression
            eigVal_LT = chebfun(@(x) ((x*x_scale+muN(x*x_scale).^2)./(x*x_scale-wN(x*x_scale).^2)+(x*x_scale-wN(x*x_scale).^2)/(x*x_scale+muN(x*x_scale).^2)) .* cos(psiN(x*x_scale)).*cos(wN(x*x_scale)) + (wN(x*x_scale)./muN(x*x_scale)-muN(x*x_scale)./wN(x*x_scale)).*sin(psiN(x*x_scale)).*sin(wN(x*x_scale))-2*gN(x*x_scale), range_LT, 'splitting', 'on');
        end
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs((-1/gammaN-gammaN) * ...
            cos(sqrt(alphaN*(1+gammaN)))+ ...
            (sqrt(alphaN*(1+gammaN)))* ...
            sin(sqrt(alphaN*(1+gammaN)))-2) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
%         eigVal_GT =  chebfun(@(x) 2*x*x_scale*(x*x_scale-wN(x*x_scale).^2).*(x*x_scale-thetaN(x*x_scale).^2)./(wN(x*x_scale).*thetaN(x*x_scale)) - x*cos(wN(x)).*cos(thetaN(x)).*(2*x^2-2*x.*wN(x).^2-2*x.*thetaN(x).^2+wN(x).^4+thetaN(x).^4)./(wN(x).*thetaN(x)) - x * sin(wN(x)).*sin(thetaN(x)).*(wN(x).^2+thetaN(x).^2).*(x-wN(x).^2).*(x-thetaN(x).^2)./(wN(x).^2.*thetaN(x).^2), ...
%             range_GT, 'splitting', 'on');
        eigVal_GT =  chebfun(@(x) ((thetaN(x*x_scale).^2-x*x_scale)./(x*x_scale-wN(x*x_scale).^2)+(x*x_scale-wN(x*x_scale).^2)./(thetaN(x*x_scale).^2-x*x_scale)) .* cos(thetaN(x*x_scale)) .* cos(wN(x*x_scale)) - (thetaN(x*x_scale).^2+wN(x*x_scale).^2)./(thetaN(x*x_scale).*wN(x*x_scale)).*sin(thetaN(x*x_scale)).*sin(wN(x*x_scale)) + 2, ...
            range_GT, 'splitting', 'on');
        %----------------------------
    case "CP"  % clamped-pinned
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT = chebfun(@(x)(muN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2).*sinh(muN(x*x_scale)).*cos(wN(x*x_scale)))./(wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2))+cosh(muN(x*x_scale)).*sin(wN(x*x_scale)), range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            eigVal_LT = chebfun( @(x) muN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)/(wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2)) .*sin(psiN(x*x_scale)).*cos(wN(x*x_scale))+cos(psiN(x*x_scale)).*sin(wN(x*x_scale)), range_LT, 'splitting', 'on');
        end
        
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(wSqrt^2 * sin(wSqrt)) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
        eigVal_GT = chebfun(@(x)(wN(x*x_scale).^2-thetaN(x*x_scale).^2).*(sin(thetaN(x*x_scale)).*cos(wN(x*x_scale))-(wN(x*x_scale).*(x*x_scale-thetaN(x*x_scale).^2).*cos(thetaN(x*x_scale)).*sin(wN(x*x_scale)))./(thetaN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2))), range_GT, 'splitting', 'on');
        %----------------------------
    case "FF"  % free-free
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT =  chebfun(@(x) -(muN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2)./(wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)) - wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)./(muN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2))).*sinh(muN(x*x_scale)).*sin(wN(x*x_scale))+2*cosh(muN(x*x_scale)).*cos(wN(x*x_scale))-2, range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            eigVal_LT = chebfun(@(x) -(muN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2)./(wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)) - wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)./(muN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2))).*sin(psiN(x*x_scale)).*sin(wN(x*x_scale)) + 2*cos(psiN(x*x_scale)).*cos(wN(x*x_scale))-2*gN(x*x_scale), range_LT, 'splitting' , 'on');
        end
        
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(alphaN^3 * gammaN/wSqrt *(2+gammaN*wSqrt*sin(wSqrt)-2*cos(wSqrt))) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
        eigVal_GT = chebfun(@(x) sin(thetaN(x*x_scale)).*sin(wN(x*x_scale)).*(-(thetaN(x*x_scale).*(x*x_scale-thetaN(x*x_scale).^2))./(2*wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2))-wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)./(2*thetaN(x*x_scale).*(x*x_scale-thetaN(x*x_scale).^2)))-cos(thetaN(x*x_scale)).*cos(wN(x*x_scale))+1, range_GT, 'splitting', 'off');
        %----------------------------
    case "PF"  % pinned-free
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT = chebfun(@(x) x*x_scale.* muN(x*x_scale).*(x*x_scale + muN(x*x_scale).^2).*sinh(muN(x*x_scale)).*cos(wN(x*x_scale)) + x*x_scale.*wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2).*cosh(muN(x*x_scale)).*sin(wN(x*x_scale)), range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            eigVal_LT = chebfun(@(x) sqrt(muN(x*x_scale)./wN(x*x_scale)) .* sin(psiN(x*x_scale)) .* cos(wN(x*x_scale)) + sqrt(wN(x*x_scale)./muN(x*x_scale)) .* (x*x_scale-wN(x*x_scale).^2)./(x*x_scale+muN(x*x_scale).^2).*cos(psiN(x*x_scale)).*sin(wN(x*x_scale)), range_LT, 'splitting', 'on');
        end
        
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(alphaN^2*gammaN*(1+gammaN)*sin(wSqrt)) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
        eigVal_GT = chebfun(@(x) (x*x_scale-wN(x*x_scale).^2)./thetaN(x*x_scale) .* cos(thetaN(x*x_scale)) .* sin(wN(x*x_scale))  - (x*x_scale-thetaN(x*x_scale).^2)./wN(x*x_scale).*cos(wN(x*x_scale)).*sin(thetaN(x*x_scale)), range_GT, 'splitting', 'on');
        %----------------------------
    case "PP"  % pinned-pinned <-check the expressions here (9/17/2018)
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT = chebfun(@(x) (-sin(wN(x*x_scale)).*(muN(x*x_scale).^2+wN(x*x_scale).^2).^2).*sinh(muN(x*x_scale)), range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            eigVal_LT = chebfun(@(x) sin(wN(x*x_scale)), range_LT, 'splitting', 'on');
        end

        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(sin(wSqrt)) < tol;
        eigVal_EQ = eq_check * alphaN;

        % lambda > alpha
        eigVal_GT = chebfun(@(x) (-sin(wN(x*x_scale)).*(wN(x*x_scale).^2-thetaN(x*x_scale).^2).^2).*sin(thetaN(x*x_scale)), range_GT, 'splitting', 'on');
        %----------------------------
    case "PR"  % pinned-roller
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT = chebfun(@(x) cos(wN(x*x_scale)).*cosh(muN(x*x_scale)), range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            eigVal_LT = chebfun(@(x) cos(wN(x*x_scale)).*cos(psiN(x*x_scale)), range_LT, 'splitting', 'on');
        end
        
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(wSqrt.^3.*cos(wSqrt)) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
        eigVal_GT = chebfun(@(x) cos(wN(x*x_scale)).*cos(thetaN(x*x_scale)), range_GT, 'splitting', 'on');
        %----------------------------
    case "CR"  % clamped-roller
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT = chebfun(@(x) sin(wN(x*x_scale)).*cosh(muN(x*x_scale))-wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2)./(muN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)).*cos(wN(x*x_scale)).*sinh(muN(x*x_scale)), range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            eigVal_LT = chebfun(@(x) sin(wN(x*x_scale)).*cos(psiN(x*x_scale))-wN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2)./(muN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)).*cos(wN(x*x_scale)).*sin(psiN(x*x_scale)), range_LT, 'splitting', 'on');
        end
        
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(-alphaN*gammaN*(sin(wSqrt)+cos(wSqrt))) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
        eigVal_GT = chebfun(@(x) sin(wN(x*x_scale)).*cos(thetaN(x*x_scale))-wN(x*x_scale).*(x*x_scale-thetaN(x*x_scale).^2)./(thetaN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)).*cos(wN(x*x_scale)).*sin(thetaN(x*x_scale)), range_GT, 'splitting', 'on');
        %----------------------------
    case "RF"  % roller-free
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT = chebfun(@(x) sin(wN(x*x_scale)).*cosh(muN(x*x_scale)) - wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)./(muN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2)).*cos(wN(x*x_scale)).*sinh(muN(x*x_scale)), range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            %             eigVal_LT = chebfun(@(x*x_scale) x*x_scale.^2.*cos(psiN(x*x_scale)).*sin(wN(x*x_scale)).*(muN(x*x_scale).^2+x*x_scale.^2).*(muN(x*x_scale).^2+wN(x*x_scale).^2)./(muN(x*x_scale).*wN(x*x_scale).^2)-x*x_scale.^2.*cos(wN(x*x_scale)).*sin(psiN(x*x_scale)).*(x*x_scale-wN(x*x_scale).^2).*(muN(x*x_scale).^2+wN(x*x_scale).^2)./(muN(x*x_scale).^2.*wN(x*x_scale)), range_LT, 'splitting', 'on');
            eigVal_LT = chebfun(@(x) sin(wN(x*x_scale)).*cos(psiN(x*x_scale)) - wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)./(muN(x*x_scale).*(x*x_scale+muN(x*x_scale).^2)).*cos(wN(x*x_scale)).*sin(psiN(x*x_scale)), range_LT, 'splitting', 'on');
            
        end
        
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(alphaN^2*(sin(wSqrt)+gammaN*wSqrt*cos(wSqrt))) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
        %         eigVal_GT = chebfun(@(x*x_scale) x*x_scale.^2.*cos(wN(x*x_scale)).*sin(thetaN(x*x_scale)).*(wN(x*x_scale).^2-thetaN(x*x_scale).^2).*(x*x_scale.^2-wN(x*x_scale).^2)./(wN(x*x_scale).*thetaN(x*x_scale).^2)-x*x_scale.^2.*cos(thetaN(x*x_scale)).*sin(wN(x*x_scale)).*(wN(x*x_scale).^2-thetaN(x*x_scale).^2).*(x*x_scale-thetaN(x*x_scale).^2)./(wN(x*x_scale).^2.*thetaN(x*x_scale)), range_GT, 'splitting', 'on');
        eigVal_GT = chebfun(@(x) sin(wN(x*x_scale)).*cos(thetaN(x*x_scale)) - wN(x*x_scale).*(x*x_scale-wN(x*x_scale).^2)./(thetaN(x*x_scale).*(x*x_scale-thetaN(x*x_scale).^2)).*cos(wN(x*x_scale)).*sin(thetaN(x*x_scale)), range_GT, 'splitting', 'on');
        %----------------------------
    case "RR"  % roller-roller
        % lambda < alpha
        if ~numStable  % this the conventional expressions
            eigVal_LT = chebfun(@(x) sin(wN(x*x_scale)).*sinh(muN(x*x_scale)), range_LT, 'splitting', 'on');
        else  % this is the stabilized expressions
            eigVal_LT = chebfun(@(x) sin(wN(x*x_scale)).*sin(psiN((x*x_scale))), range_LT, 'splitting' ,'on');
        end
        
        % check if the case alpha = lambda, within tol, exists
        eq_check = abs(alphaN^2*(1+gammaN)*sin(wSqrt)) < tol;
        eigVal_EQ = eq_check * alphaN;
        
        % lambda > alpha
        eigVal_GT = chebfun(@(x) sin(wN(x*x_scale)).*sin(thetaN(x*x_scale)), range_GT, 'splitting', 'on');
end
% prepare the output
cheb_rep = {eigVal_LT, eigVal_EQ(eigVal_EQ > 0), eigVal_GT};

end