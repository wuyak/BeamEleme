% function cheb_rep = getTimoshenkoModes(sys, BC, lambda_cell)
% A function for finding the eigenmodes for the Timoshenko beam.
% it returns a chebfun represenation that can be qeried for roots, plotted,
% differentiated, or integrated.
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
% lambda_cell: a 3x1 cell area which contains the eigenvalues in the three
%              different regimes lambda<alpha, lambda=alpha, lambda>alpha,
%              respecitvely. These roots, e.g., can be obtained using
%              roots() with the output of getTimoshenkoEigenvalues.m.
%
% Output:
% chep_rep: a nested 3x1 cell array of chebfunction representations. Each
% cell contains the output that corresponds to lambda<alpha, lambda=alpha,
% and lambda>alpha, respectively. Each one of those cells further contains
% another cell with 4 entries that contains {U, phi, U', phi'}. In summary,
% the output is:
% cheb_rep = {{eigmode_U_LT, eigmode_phi_LT, eigmode_U_LT_prime, eigmode_phi_LT_prime}, ...
%     {eigmode_U_EQ, eigmode_phi_EQ, eigmode_U_EQ_prime, eigmode_phi_EQ_prime}, ...
%     {eigmode_U_GT, eigmode_phi_GT, eigmode_U_GT_prime, eigmode_phi_GT_prime}};
%
function cheb_rep = getTimoshenkoModes(sys, BC, lambda_cell)

% Add chebfun to path if needed
if ~exist('chebfun', 'file')
    chebfun_path = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'chebfun');
    if exist(chebfun_path, 'dir')
        addpath(genpath(chebfun_path));
    else
        error('chebfun not found. Please install from https://www.chebfun.org/ or check path: %s', chebfun_path);
    end
end

% adjust for lambda scaling, if any
if isfield(sys, 'lambda_scale')
    lambda_cell =  cellfun(@(piece) sys.lambda_scale * piece, lambda_cell, 'UniformOutput', false);
end

% extract the lambdas in the different regimes
lambda_LT = lambda_cell{1};
lambda_LT = lambda_LT.';  % turn into a row vector
%
lambda_EQ = lambda_cell{2};
lambda_EQ = lambda_EQ.';  % turn into a row vector
%
lambda_GT = lambda_cell{3};
lambda_GT = lambda_GT.';  % turn into a row vector

% make the BC string upper case
BC = upper(BC);

% extract needed parameters
alphaN = sys.par.alpha;
gammaN = sys.par.gamma;
% betaN = sys.par.beta;
% maxMode = sys.maxMode;
numStable = sys.numStable;
% wSqrt = sqrt(alphaN*(1+gammaN));

% get the needed useful params
% get the needed useful params
deltaN = @(x) 1 - 4*gammaN/(1+gammaN).^2 .* (1-alphaN./x);
muN = @(x) sqrt(0.5*x.*(1+gammaN).*(sqrt(deltaN(x))-1));
wN = @(x) sqrt(0.5*x*(1+gammaN).*(sqrt(deltaN(x))+1));
thetaN = @(x) sqrt(0.5*x*(1+gammaN).*(1-sqrt(deltaN(x))));
% psiN = @(x) atan2(1-exp(-2*muN(x)), 1+exp(-2*muN(x)));
% gN = @(x) sqrt(2)*exp(-muN(x))./sqrt(1+exp(-4*muN(x)));

% set the nondiemnsional beam span
z_range = [0, 1];  % non-diemensional beam span

%% pick the correct expressions based on the BC and the lambda case
% initiate the flag for alternative stable expressions
flag_alt1_stab = false;
switch BC
    case "CC"    % fixed-fixed
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = ones(1, length(lambda_LT));
                A2 = (-mu1Vec.*(lambda_LT-w1Vec.^2).*sinh(mu1Vec)-w1Vec.*(lambda_LT+mu1Vec.^2).*sin(w1Vec))./...
                    ((mu1Vec.*(lambda_LT-w1Vec.^2)).*cosh(mu1Vec) - mu1Vec.*(lambda_LT-w1Vec.^2).*cos(w1Vec));
                A3 = w1Vec.*(lambda_LT + mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                A4 = -A2;
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                P = 2*(-exp(-mu1Vec) + w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2)) .* sin(w1Vec)+cos(w1Vec))./...
                    (1+exp(-2*mu1Vec)-2*exp(-mu1Vec).*cos(w1Vec)) ;
                A3_stab = w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                R = -1;
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = ones(1, length(lambda_EQ));
            A2 = w2Vec.*sin(w2Vec)./((cos(w2Vec)-1).*(lambda_EQ-w2Vec.^2));
            A3 = w2Vec./(lambda_EQ-w2Vec.^2);
            A4 = -A2;
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = (theta3Vec.*(lambda_GT - w3Vec.^2).*sin(theta3Vec)-w3Vec.*(lambda_GT-theta3Vec.^2).*sin(w3Vec))./ ...
                ((theta3Vec.*(lambda_GT-w3Vec.^2)).*cos(w3Vec) -(theta3Vec.*(lambda_GT-w3Vec.^2)).*cos(theta3Vec) );
            A3 = -w3Vec.*(lambda_GT-theta3Vec.^2)./(theta3Vec.*(lambda_GT-w3Vec.^2));
            A4 = -A2;
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "CF"    % fixed-free (cantilever beam)
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = ones(1, length(lambda_LT));
                A2 = -(w1Vec.*cos(w1Vec).*(mu1Vec.^2+lambda_LT)-w1Vec.*cosh(mu1Vec).*(lambda_LT-w1Vec.^2))./...
                    (mu1Vec.*sin(w1Vec).*(lambda_LT-w1Vec.^2)-w1Vec.*sinh(mu1Vec).*(lambda_LT-w1Vec.^2));
                A3 = w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                A4 = -A2;
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                P = 2*(-w1Vec.*exp(-mu1Vec).*(lambda_LT-w1Vec.^2)+w1Vec.*(lambda_LT+mu1Vec.^2).*cos(w1Vec)-mu1Vec.*(lambda_LT-w1Vec.^2).*sin(w1Vec))./...
                    ((lambda_LT-w1Vec.^2).*(w1Vec.*(exp(-2*mu1Vec)-1)+2*mu1Vec.*exp(-mu1Vec) .*sin(w1Vec))   );
                %                 P = (2*w1Vec.^3.*exp(-mu1Vec)+(2*w1Vec.*mu1Vec.^2+2*lambda_LT.*w1Vec).*cos(w1Vec)+(2*mu1Vec.*w1Vec.^2-2*lambda_LT.*mu1Vec).*sin(w1Vec)-2*lambda_LT.*w1Vec.*exp(-mu1Vec))./ ...
                %                     ((2*mu1Vec.*exp(-mu1Vec).*(lambda_LT-w1Vec.^2)).*sin(w1Vec)+w1Vec.*(exp(-1*mu1Vec)-1).*(lambda_LT-w1Vec.^2));
                A3_stab = w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                R = -1;
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = ones(1, length(lambda_EQ));
            A2 = (w2Vec.*(lambda_EQ.*cos(w2Vec)-lambda_EQ+w2Vec.^2))./...
                ((alphaN * w2Vec-lambda_EQ.*sin(w2Vec)).*(lambda_EQ-w2Vec.^2));
            A3 = w2Vec./(lambda_EQ-w2Vec.^2);
            A4 = -A2;
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = -(((lambda_GT-theta3Vec.^2)./(lambda_GT-w3Vec.^2)).*cos(w3Vec)-cos(theta3Vec))./...
                (sin(theta3Vec)- (theta3Vec./w3Vec) .* sin(w3Vec));
            A3 = -w3Vec.*(lambda_GT-theta3Vec.^2)./(theta3Vec.*(lambda_GT-w3Vec.^2));
            A4 = -A2;
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "CP"  % clamped-pinned
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = ones(1, length(lambda_LT));
                A2 = -(sinh(mu1Vec) + w1Vec.*(lambda_LT+mu1Vec.^2) ./ (mu1Vec.*(lambda_LT-w1Vec.^2)) .*sin(w1Vec))./...
                    (cosh(mu1Vec)-cos(w1Vec));
                A3 = w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                A4 = -A2;
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                P = 2*(-mu1Vec.*exp(-mu1Vec).*(lambda_LT-w1Vec.^2)+w1Vec.*(lambda_LT+mu1Vec.^2).*sin(w1Vec)+mu1Vec.*(lambda_LT-w1Vec.^2).*cos(w1Vec))./...
                    (mu1Vec.*(exp(-2*mu1Vec)+1).*(lambda_LT-w1Vec.^2)-2*mu1Vec.*exp(-mu1Vec).*(lambda_LT-w1Vec.^2).*cos(w1Vec));
                A3_stab = w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                R = -1;
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = ones(1, length(lambda_EQ));
            A2 = w2Vec.*sin(w2Vec)./((lambda_EQ-w2Vec.^2).*(cos(w2Vec)-1));
            A3 = w2Vec./(lambda_EQ-w2Vec.^2);
            A4 = -A2;
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = (sin(theta3Vec)-w3Vec.*(lambda_GT-theta3Vec.^2)./(theta3Vec.*(lambda_GT-w3Vec.^2)).*sin(w3Vec))./...
                (cos(w3Vec)-cos(theta3Vec));
            A3 = -w3Vec.*(lambda_GT-theta3Vec.^2)./(theta3Vec.*(lambda_GT-w3Vec.^2));
            A4 = -A2;
            %
            A_GT = [A1; A2; A3; A4];
        end
        
        %----------------------------
    case "FF"  % free-free
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = ones(1, length(lambda_LT));
                A2 = -(cosh(mu1Vec)-cos(w1Vec))./...
                    (sinh(mu1Vec)-((mu1Vec.*(lambda_LT+mu1Vec.^2))./(w1Vec.*(lambda_LT-w1Vec.^2))).*sin(w1Vec));
                A3 = w1Vec./mu1Vec;
                A4 = (cosh(mu1Vec)-cos(w1Vec))./...
                    ((lambda_LT-w1Vec.^2)./(lambda_LT+mu1Vec.^2)).*(sinh(mu1Vec)-(mu1Vec./w1Vec).*sin(w1Vec));
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                P = -2*(w1Vec.*exp(-mu1Vec).*(lambda_LT-w1Vec.^2)+mu1Vec.*(lambda_LT+mu1Vec.^2).*sin(w1Vec)-w1Vec.*(lambda_LT-w1Vec.^2).*cos(w1Vec))./...
                    (2*mu1Vec.*exp(-mu1Vec).*(lambda_LT+mu1Vec.^2).*sin(w1Vec)+w1Vec.*(exp(-2*mu1Vec)-1).*(lambda_LT-w1Vec.^2));
                A3_stab = w1Vec./mu1Vec;
                R = -(lambda_LT+mu1Vec.^2)./(lambda_LT-w1Vec.^2);
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = ones(1, length(lambda_EQ));
            A2 = -(w2Vec.*(cos(w2Vec)-1).*(lambda_EQ-w2Vec.^2))./...
                (alphaN*(lambda_EQ.*sin(w2Vec)-lambda_EQ.*w2Vec+w2Vec.^3));
            A3 = w2Vec/lambda_EQ;
            A4 = w2Vec.*(cos(w2Vec)-1)./...
                (lambda_EQ.*sin(w2Vec)-lambda_EQ.*w2Vec+w2Vec.^3);
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = -(cos(w3Vec)-cos(theta3Vec))./...
                (sin(theta3Vec)-(theta3Vec.*(lambda_GT-theta3Vec.^2)./(w3Vec.*(lambda_GT-w3Vec.^2))).*sin(w3Vec));
            A3 = -w3Vec./theta3Vec;
            A4 = (cos(w3Vec)-cos(theta3Vec))./...
                (((lambda_GT-w3Vec.^2)./(lambda_GT-theta3Vec.^2)).*sin(theta3Vec)-(theta3Vec./w3Vec).*sin(w3Vec));
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "PF"  % pinned-free
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = mu1Vec.*cos(w1Vec)./(w1Vec.*cosh(mu1Vec));
                A2 = zeros(1, length(lambda_LT));
                A3 = ones(1, length(lambda_LT));
                A4 = zeros(1, length(lambda_LT));
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                % indicate alternative procedure for stable expressions exist
                flag_alt1_stab = true;
                % write down the alternative stable expressions
                eigmode_U_LT = chebfun(@(z) (mu1Vec.*cos(w1Vec)./w1Vec).*(exp(-mu1Vec.*(1-z))-exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec)) + sin(z*w1Vec), z_range);
                eigmode_phi_LT = chebfun(@(z) (lambda_LT+mu1Vec.^2).*cos(w1Vec)./w1Vec .* (exp(-mu1Vec.*(1-z))+exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec)) - (lambda_LT-w1Vec.^2)./w1Vec .* cos(z*w1Vec), z_range);
                % derivatives
                eigmode_U_LT_prime = chebfun(@(z) (mu1Vec./w1Vec.*cos(w1Vec)).*mu1Vec.*(exp(-mu1Vec.*(1-z))+exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec)) + w1Vec .* cos(z*w1Vec), z_range);
                eigmode_phi_LT_prime = chebfun(@(z) (lambda_LT+mu1Vec.^2)./w1Vec .* cos(w1Vec) .* mu1Vec .* (exp(-mu1Vec.*(1-z)) - exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec)) + (lambda_LT-w1Vec.^2) .* sin(z*w1Vec), z_range);
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = ones(1, length(lambda_EQ));
            A2 = zeros(1, length(lambda_EQ));
            A3 = w2Vec./(lambda_EQ.*cos(w2Vec));
            A4 = zeros(1, length(lambda_EQ));
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = zeros(1, length(lambda_GT));
            A3 = -w3Vec.*cos(theta3Vec)./(theta3Vec.*cos(w3Vec));
            A4 = zeros(1, length(lambda_GT));
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "PP"  % pinned-pinned
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = -sin(w1Vec)./sinh(mu1Vec);
                A2 = zeros(1, length(lambda_LT));
                A3 = ones(1, length(lambda_LT));
                A4 = zeros(1, length(lambda_LT));
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                % indicate alternative procedure for stable expressions exist
                flag_alt1_stab = true;
                % write down the alternative stable expressions
                eigmode_U_LT = chebfun(@(z) sin(w1Vec.*z), z_range);
                eigmode_phi_LT = chebfun(@(z) -(lambda_LT-w1Vec.^2)./w1Vec.*cos(w1Vec.*z), z_range);
                % derivatives
                eigmode_U_LT_prime = chebfun(@(z) w1Vec.*cos(w1Vec.*z), z_range);
                eigmode_phi_LT_prime = chebfun(@(z) (lambda_LT-w1Vec.^2) .* sin(w1Vec.*z), z_range);
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            % indicate alternative procedure for stable expressions exist
            flag_alt2_stab = true;

            w2Vec = wN(alphaN);
            %
            eigmode_U_EQ = chebfun(@(z) sin(w2Vec.*z), z_range);
            eigmode_phi_EQ = chebfun(@(z) -(lambda_EQ-w2Vec.^2)./w2Vec.*cos(w2Vec.*z), z_range);
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = zeros(1, length(lambda_GT));
            A3 = -sin(theta3Vec)./sin(w3Vec);
            A4 = zeros(1, length(lambda_GT));
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "PR"  % pinned-roller
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = cos(w1Vec).*(mu1Vec.*(lambda_LT-w1Vec.^2))./(w1Vec.*cosh(mu1Vec).*(lambda_LT+mu1Vec.^2));
                A2 = zeros(1, length(lambda_LT));
                A3 = ones(1, length(lambda_LT));
                A4 = zeros(1, length(lambda_LT));
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                % indicate alternative procedure for stable expressions exist
                flag_alt1_stab = true;
                % write down the alternative stable expressions
                eigmode_U_LT = chebfun(@(z) sin(w1Vec.*z), z_range);
                eigmode_phi_LT = chebfun(@(z) -(lambda_LT-w1Vec.^2)./w1Vec.*cos(w1Vec.*z), z_range);
                % derivatives
                eigmode_U_LT_prime = chebfun(@(z) w1Vec .* cos(w1Vec.*z), z_range);
                eigmode_phi_LT_prime = chebfun(@(z) (lambda_LT-w1Vec.^2) .* sin(w1Vec.*z), z_range);
                
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = ones(1, length(lambda_EQ));
            A2 = zeros(1, length(lambda_EQ));
            A3 = w2Vec./((lambda_EQ-w2Vec.^2).*cos(w2Vec));
            A4 = zeros(1, length(lambda_EQ));
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = zeros(1, length(lambda_GT));
            A3 = -cos(theta3Vec).*(w3Vec.*(lambda_GT-theta3Vec.^2))./ ...
                (theta3Vec.*(lambda_GT-w3Vec.^2).*cos(w3Vec));
            A4 = zeros(1, length(lambda_GT));
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "CR"  % clamped-roller
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = ones(1, length(lambda_LT));
                A2 = (cosh(mu1Vec)-cos(w1Vec))./...
                    ((mu1Vec.*(lambda_LT-w1Vec.^2)./(w1Vec.*(lambda_LT+mu1Vec.^2))).*sin(w1Vec)-sinh(mu1Vec));
                A3 = w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                A4 = -A2;
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                P = 2*(w1Vec.*(lambda_LT+mu1Vec.^2).*cos(w1Vec)-mu1Vec.*(lambda_LT-w1Vec.^2).*sin(w1Vec)-w1Vec.*(lambda_LT+mu1Vec.^2).*exp(-mu1Vec))./...
                    (-w1Vec.*(lambda_LT+mu1Vec.^2)+2*mu1Vec.*(lambda_LT-w1Vec.^2).*exp(-mu1Vec).*sin(w1Vec)+w1Vec.*(lambda_LT+mu1Vec.^2).*exp(-2*mu1Vec));
                A3_stab = w1Vec.*(lambda_LT+mu1Vec.^2)./(mu1Vec.*(lambda_LT-w1Vec.^2));
                R = -1;
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = ones(1, length(lambda_EQ));
            A2 = -w2Vec.*(cos(w2Vec)-1)./((lambda_EQ-w2Vec.^2).*sin(w2Vec)-alphaN*w2Vec);
            A3 = w2Vec./(lambda_EQ-w2Vec.^2);
            A4 = w2Vec.*(1-cos(w2Vec))./((lambda_EQ-w2Vec.^2).*sin(w2Vec)+alphaN*w2Vec);
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = ones(1, length(lambda_GT));
            A2 = (cos(w3Vec)-cos(theta3Vec)) ./ ...
                ((theta3Vec.*(lambda_GT-w3Vec.^2)./(w3Vec.*(lambda_GT-theta3Vec.^2))).*sin(w3Vec)-sin(theta3Vec));
            A3 = -(w3Vec.*(lambda_GT-theta3Vec.^2))./(theta3Vec.*(lambda_GT-w3Vec.^2));
            A4 = -A2;
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "RF"  % roller-free
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = zeros(1, length(lambda_LT));
                A2 = -(lambda_LT-w1Vec.^2).*cos(w1Vec)./((lambda_LT+mu1Vec.^2).*cosh(mu1Vec));
                A3 = zeros(1, length(lambda_LT));
                A4 = ones(1, length(lambda_LT));
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                % indicate alternative procedure for stable expressions exist
                flag_alt1_stab = true;
                % write down the alternative stable expressions
                eigmode_U_LT = chebfun(@(z) -((lambda_LT-w1Vec.^2).*cos(w1Vec)./(lambda_LT+mu1Vec.^2)).*(exp(-mu1Vec.*(1-z))+exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec))+cos(w1Vec.*z), z_range);
                eigmode_phi_LT = chebfun(@(z) -((lambda_LT-w1Vec.^2).*cos(w1Vec)./mu1Vec).*(exp(-mu1Vec.*(1-z))-exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec))+(lambda_LT-w1Vec.^2)./w1Vec .* sin(w1Vec.*z), z_range);
                % derivatives
                eigmode_U_LT_prime = chebfun(@(z) -((lambda_LT-w1Vec.^2).*cos(w1Vec)./(lambda_LT+mu1Vec.^2)).*mu1Vec.*(exp(-mu1Vec.*(1-z))-exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec))-w1Vec.*sin(w1Vec.*z), z_range);
                eigmode_phi_LT_prime = chebfun(@(z) -((lambda_LT-w1Vec.^2).*cos(w1Vec)).*(exp(-mu1Vec.*(1-z)) + exp(-mu1Vec.*(1+z)))./(1+exp(-2*mu1Vec))+(lambda_LT-w1Vec.^2) .* cos(w1Vec.*z), z_range);
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = zeros(1, length(lambda_EQ));
            A2 = ones(1, length(lambda_EQ));
            A3 = zeros(1, length(lambda_EQ));
            A4 = -alphaN./((lambda_EQ-w2Vec.^2).*cos(w2Vec));
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = zeros(1, length(lambda_GT));
            A2 = ones(1, length(lambda_GT));
            A3 = zeros(1, length(lambda_GT));
            A4 = -(lambda_GT-theta3Vec.^2).*cos(theta3Vec)./((lambda_GT-w3Vec.^2).*cos(w3Vec));
            %
            A_GT = [A1; A2; A3; A4];
        end
        %----------------------------
    case "RR"  % roller-roller
        % lambda < alpha
        % get the coefficients obtained from the eigenvalue problem. if the
        % stabilized epxressions are requested, get P, A3, and R instead.
        if ~isempty(lambda_LT)
            mu1Vec = muN(lambda_LT);
            w1Vec = wN(lambda_LT);
            if ~numStable
                A1 = zeros(1, length(lambda_LT));
                A2 = -mu1Vec.*(lambda_LT-w1Vec.^2).*sin(w1Vec)./(w1Vec.*(lambda_LT+mu1Vec.^2).*sinh(mu1Vec));
                A3 = zeros(1, length(lambda_LT));
                A4 = ones(1, length(lambda_LT));
                
                % put the terms into one matrix
                A_LT = [A1; A2; A3; A4];
            else
                % indicate alternative procedure for stable expressions exist
                flag_alt1_stab = true;
                % write down the alternative stable expressions
                eigmode_U_LT = chebfun(@(z) cos(w1Vec.*z), z_range);
                eigmode_phi_LT = chebfun(@(z) (lambda_LT-w1Vec.^2)./w1Vec.*sin(w1Vec.*z), z_range);
                % derivatives
                eigmode_U_LT_prime = chebfun(@(z) -w1Vec.* sin(w1Vec.*z), z_range);
                eigmode_phi_LT_prime = chebfun(@(z) (lambda_LT-w1Vec.^2) .* cos(w1Vec.*z), z_range);
            end
        end
        % lambda = alpha
        if ~isempty(lambda_EQ)
            w2Vec = wN(alphaN);
            A1 = zeros(1, length(lambda_EQ));
            A2 = ones(1, length(lambda_EQ));
            A3 = zeros(1, length(lambda_EQ));
            A4 = -alphaN./((lambda_EQ-w2Vec.^2).*cos(w2Vec));
            %
            A_EQ = [A1; A2; A3; A4];
        end
        % lambda > alpha
        if ~isempty(lambda_GT)
            w3Vec = wN(lambda_GT);
            theta3Vec = thetaN(lambda_GT);
            %
            A1 = zeros(1, length(lambda_GT));
            A2 = ones(1, length(lambda_GT));
            A3 = zeros(1, length(lambda_GT));
            A4 = -alphaN*w3Vec./((lambda_GT-w3Vec.^2).*sin(w3Vec));
            %
            A_GT = [A1; A2; A3; A4];
        end
end

%% prepare the output
% modes for lambda < alpha
if ~isempty(lambda_LT)
    if ~numStable  % these are the conventional eigenmodes
        eigmode_U_LT = chebfun(@(z) ...
            A_LT(1, :) .* sinh(z*mu1Vec) + ...
            A_LT(2, :) .* cosh(z*mu1Vec) + ...
            A_LT(3, :) .* sin(z*w1Vec)   + ...
            A_LT(4, :) .* cos(z*w1Vec), ...
            z_range, 'splitting', 'on');
        eigmode_phi_LT = chebfun(@(z) ...
            A_LT(1, :) .* (lambda_LT+mu1Vec.^2)./mu1Vec .* cosh(z*mu1Vec) +...
            A_LT(2, :) .*  (lambda_LT+mu1Vec.^2)./mu1Vec .* sinh(z*mu1Vec) + ...
            A_LT(3, :) .* -(lambda_LT-w1Vec.^2)./w1Vec .* cos(z*w1Vec)   + ...
            A_LT(4, :) .*  (lambda_LT-w1Vec.^2)./w1Vec  .* sin(z*w1Vec), ...
            z_range, 'splitting', 'on');
        % derivatives
        eigmode_U_LT_prime = chebfun(@(z) ...
            A_LT(1, :) .* mu1Vec.*cosh(z*mu1Vec) + ...
            A_LT(2, :) .* mu1Vec.*sinh(z*mu1Vec) + ...
            A_LT(3, :) .* w1Vec .* cos(z*w1Vec)   + ...
            A_LT(4, :) .* w1Vec .* -sin(z*w1Vec), ...
            z_range, 'splitting', 'on');
        eigmode_phi_LT_prime = chebfun(@(z) ...
            A_LT(1, :) .*  (lambda_LT+mu1Vec.^2) .* sinh(z*mu1Vec) +...
            A_LT(2, :) .*  (lambda_LT+mu1Vec.^2) .* cosh(z*mu1Vec) + ...
            A_LT(3, :) .* (lambda_LT-w1Vec.^2) .* sin(z*w1Vec)   + ...
            A_LT(4, :) .*  (lambda_LT-w1Vec.^2) .* cos(z*w1Vec), ...
            z_range, 'splitting', 'on');
        % these are the stabilized modes
        % if alternative epxressions were used (flag_alt_stab=1), skip this step
    elseif ~flag_alt1_stab
        eigmode_U_LT = chebfun(@(z) -exp(-z.*mu1Vec) - P/2.*(exp(-(1-z).*mu1Vec) + exp(-(1+z).*mu1Vec)) + A3_stab .* sin(z.*w1Vec) - R .* (1+P .* exp(-mu1Vec)).*cos(z.*w1Vec), ...
            z_range, 'splitting', 'on');
        eigmode_phi_LT = chebfun(@(z) (lambda_LT+mu1Vec.^2)./mu1Vec .*(exp(-z.*mu1Vec) - P/2 .*(exp(-(1-z).*mu1Vec) - exp(-(1+z).*mu1Vec)))-A3_stab .* ((lambda_LT-w1Vec.^2)./w1Vec).*cos(z.*w1Vec)-R.*(1+P.*exp(-mu1Vec)) .* ((lambda_LT-w1Vec.^2)./w1Vec).* sin(z.*w1Vec), ...
            z_range, 'splitting', 'on');
        % derivatives
        eigmode_U_LT_prime = chebfun(@(z) mu1Vec.*exp(-z.*mu1Vec) - P/2.* mu1Vec.*(exp(-(1-z).*mu1Vec) - exp(-(1+z).*mu1Vec)) + A3_stab .* w1Vec.*cos(z.*w1Vec) + R .* (1+P .* exp(-mu1Vec)).*w1Vec.*sin(z.*w1Vec), ...
            z_range, 'splitting', 'on');
        eigmode_phi_LT_prime = chebfun(@(z) (lambda_LT+mu1Vec.^2) .*(-exp(-z.*mu1Vec) - P/2 .* (exp(-(1-z).*mu1Vec) + exp(-(1+z).*mu1Vec))) + A3_stab .* (lambda_LT-w1Vec.^2).*sin(z.*w1Vec)-R.*(1+P.*exp(-mu1Vec)) .*(lambda_LT-w1Vec.^2) .* cos(z.*w1Vec), ...
            z_range, 'splitting', 'on');
    end
else
    eigmode_U_LT = [];
    eigmode_phi_LT = [];
    eigmode_U_LT_prime = [];
    eigmode_phi_LT_prime = [];
end

% modes for lambda = alpha
if ~isempty(lambda_EQ)
    if ~flag_alt2_stab
        eigmode_U_EQ = chebfun(@(z) ...
            A_EQ(2, :) + ...
            A_EQ(3, :) .* sin(z*w2Vec) + ...
            A_EQ(4, :) .* cos(z*w2Vec), ...
            z_range, 'splitting', 'on');
        eigmode_phi_EQ = chebfun(@(z) ...
            A_EQ(1, :) +...
            A_EQ(2, :) .*  alphaN*z + ...
            A_EQ(3, :) .* -(lambda_EQ-w2Vec.^2)./w2Vec.*cos(z*w2Vec) + ...
            A_EQ(4, :) .* (lambda_EQ-w2Vec.^2)./w2Vec .* sin(z*w2Vec), ...
            z_range, 'splitting', 'on');
        % derivatives
        eigmode_U_EQ_prime = chebfun(@(z) ...
            A_EQ(3, :) .* w2Vec .* cos(z*w2Vec) + ...
            A_EQ(4, :) * w2Vec.* -sin(z*w2Vec), ...
            z_range, 'splitting', 'on');
        eigmode_phi_EQ_prime = chebfun(@(z) ...
            A_EQ(2, :) .*  alphaN + ...
            A_EQ(3, :) .* (lambda_EQ-w2Vec.^2).*sin(z*w2Vec) + ...
            A_EQ(4, :) .* (lambda_EQ-w2Vec.^2) .* cos(z*w2Vec), ...
            z_range, 'splitting', 'on');
    end
else
    eigmode_U_EQ = [];
    eigmode_phi_EQ = [];
    eigmode_U_EQ_prime = [];
    eigmode_phi_EQ_prime= [];
end

% modes for lambda > alpha
if ~isempty(lambda_GT)
    eigmode_U_GT =  chebfun(@(z) ...
        A_GT(1, :) .* sin(z * theta3Vec) + ...
        A_GT(2, :) .* cos(z * theta3Vec) + ...
        A_GT(3, :) .* sin(z * w3Vec) + ...
        A_GT(4, :) .* cos(z * w3Vec), ...
        z_range, 'splitting', 'off');
    eigmode_phi_GT = chebfun(@(z) ...
        A_GT(1, :) .* -(lambda_GT-theta3Vec.^2)./theta3Vec .* cos(z * theta3Vec) + ...
        A_GT(2, :) .*  (lambda_GT-theta3Vec.^2)./theta3Vec .* sin(z * theta3Vec) + ...
        A_GT(3, :) .* -(lambda_GT-w3Vec.^2)./w3Vec .* cos(z * w3Vec) + ...
        A_GT(4, :) .*  (lambda_GT-w3Vec.^2)./w3Vec .* sin(z * w3Vec), ...
        z_range, 'splitting', 'on');
    % derivative
    eigmode_U_GT_prime =  chebfun(@(z) ...
        A_GT(1, :) .* theta3Vec .* cos(z * theta3Vec) + ...
        A_GT(2, :) .* theta3Vec .* -sin(z * theta3Vec) + ...
        A_GT(3, :) .* w3Vec .* cos(z * w3Vec) + ...
        A_GT(4, :) .* w3Vec .* -sin(z * w3Vec), ...
        z_range, 'splitting', 'off');
    
    eigmode_phi_GT_prime = chebfun(@(z) ...
        A_GT(1, :) .* (lambda_GT-theta3Vec.^2) .* sin(z * theta3Vec) + ...
        A_GT(2, :) .*  (lambda_GT-theta3Vec.^2) .* cos(z * theta3Vec) + ...
        A_GT(3, :) .* (lambda_GT-w3Vec.^2) .* sin(z * w3Vec) + ...
        A_GT(4, :) .*  (lambda_GT-w3Vec.^2) .* cos(z * w3Vec), ...
        z_range, 'splitting', 'on');
else
    eigmode_U_GT = [];
    eigmode_phi_GT = [];
    eigmode_U_GT_prime = [];
    eigmode_phi_GT_prime = [];
end

cheb_rep = {{eigmode_U_LT, eigmode_phi_LT, eigmode_U_LT_prime, eigmode_phi_LT_prime}, ...
    {eigmode_U_EQ, eigmode_phi_EQ, eigmode_U_EQ_prime, eigmode_phi_EQ_prime}, ...
    {eigmode_U_GT, eigmode_phi_GT, eigmode_U_GT_prime, eigmode_phi_GT_prime}};

% normalize the modes
cheb_rep = normalize_TB_modes(sys, cheb_rep);
end

%% normalize the mode shapes
function normalized_modes = normalize_TB_modes(sys, modes)
% see if a normalization was requested
if ~isfield(sys, 'norm') || strcmp(sys.norm, 'unscaled')
    % fprintf('No normalization is requested, or incorrect or no normalization option is set. The modes are not normalized.\n');
    normalized_modes = modes;
else
    requested_norm = sys.norm;  % get the requested normalization
    % preallocate the scaled modes
    normalized_modes = cell(1, 3);
    switch requested_norm
        case "mass"  % mass orthonormalization
            for case_no = 1:3
                if ~isempty(modes{case_no}{1})
                    scale_factor = 1./sum(modes{case_no}{1}.^2+1/sys.par.alpha*modes{case_no}{2}.^2);
                    normalized_modes{case_no} = cellfun(@(x) scale_factor .* x, modes{case_no}, 'UniformOutput', false);
                else
                    normalized_modes{case_no} = modes{case_no};
                end
                
            end
    end
end
end