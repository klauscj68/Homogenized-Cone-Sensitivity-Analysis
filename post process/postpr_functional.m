function [functionals] = postpr_functional(tact,trec,taxis,...
                                           J_tot,J_drop,...
                                           mass_Evol,...
                                           datamat,pointer,...
                                           flag_Sbl)
%Given trial data set, compute the value of all associated functionals
%   tact, trec are 2 x 1 and give time ranges for activation and recovery
%   J_tot,J_drop are outputs of 'Msg' in simulation
%   mass_Evol is total number of E* across volume
%   flag_Sbl tells the script whether you need to for cG_d, Ca_d, J_d being
%    in admissible ranges (flag = false) or if you are not worried about it
%    (flag = true). If this argument is not passed, code assumes flag_Sbl =
%     false.
%   NOTE: It is expected that section L2 FIT ERROR is an exact copy of the
%         matching fitness.m value SAVE that it returns NaN if passed
%         parameters do not satisfy wellposed_ss criteria
%% Initialize
functionals = zeros(16,1);

%% Iact
% Find time points used in activation
logind = logical(...
                (tact(1) <= taxis).*(taxis <= tact(2))...
                 );
tpts = taxis(logind);
Ipts = J_tot(1) - J_tot(logind);
                             

% Find best fitting (1/2)A_I t^2 quadratic
A = bestfit(tpts,Ipts,'quad');

% Store
functionals(1) = A;

%% Ipeak
peak = max(J_drop);

if peak < 0
    peak = 0;
elseif peak > 100
    peak = 100;
end

% Store
functionals(2) = peak;

%% Irec
logind = logical(...
                (trec(1) <= taxis).*(taxis <= trec(2))...
                 );
tpts = taxis(logind);
Ipts = J_tot(1) - J_tot(logind);

% If there was overshoot then rising response, restrict to part of curve
%   strictly before the minimum and have the minimum be the 0 value. 
%   Note: Is a boundary case where minimum is at end but curve still
%         falling where this approach would return an adjusted exp
%         constant. Expect this would only occur rarely
[bottom,id] = min(Ipts);
if (bottom <= 0)
    if id == length(Ipts)
        warning('WIpts: J min was neg and at end of simulation time');
    end
    tpts = tpts(1:id-1);
    Ipts = Ipts(1:id-1);
    Ipts = Ipts - bottom;
end

try
    % Compute the c*exp(-alpha*t) best fit where min of Ipts is set to 0
    [alpha,~] = bestfit(tpts,Ipts,'exp');

    % Store
    functionals(3) = alpha;
catch ME
    warning('WIrec: There was an error with the fit');
    disp(ME.identifier);
    functionals(3) = NaN;
end


%% Estact
%   Find the points in the activation phase
logind = logical(...
                (tact(1) <= taxis).*(taxis <= tact(2))...
                 );
tpts = taxis(logind);
Epts = mass_Evol(logind);

%   Compute the quadratic .5*A*t^2 best fit
A = bestfit(tpts,Epts,'quad');

% Store
functionals(4) = A;

%% Estpeak 
Epeak = max(mass_Evol);

if Epeak > 2500
    Epeak = 2500;
end

% Store
functionals(5) = Epeak;
%% Estrec
logind = logical(...
                (trec(1) <= taxis).*(taxis <= trec(2))...
                 );
tpts = taxis(logind);
Epts = mass_Evol(logind);

%   Compute the c*exp(-alpha*t) best fit
try
    [alpha,~] = bestfit(tpts,Epts,'exp');

    % Store
    functionals(6) = alpha;
catch ME
    warning('WErec: There was an error with the fit');
    disp(ME.identifier);
    functionals(6) = NaN;
end

%% L2 Fit Error
% If Sbl trial don't worry about checking for solvability as guaranteed by
%  ranges in mc_sample_Sbl
try
    switch flag_Sbl
        case true
            flag_proceed = true;
        case false
            % Check if passed parameters are admissible by wellposed_ss
            [flag_proceed,~] = wellposed_ss(datamat,pointer);
        otherwise
            error('flag_Sbl must be true or false')
    end
catch
    % Check if passed parameters are admissible by wellposed_ss
    [flag_proceed,~] = wellposed_ss(datamat,pointer);
end

switch flag_proceed
    case true
    % Reason for these value choices is the output of flashimg
    %  Interval is range on coordinates of flashimg
    t_interval = [0.0177 0.35];
    %  Experimental trace placed flash at 5ms
    tshift = .005; 

    % Broke the image import up into three subintervals
    %  bc fitting routine was more robust this way
    p1 = @(x) (-2.509e5)*x.^3+(1.563e4)*x.^2+206.5*x-5.33;
    p2 = @(x) (141.3)*x.^2-194.9*x+21.74;
    p3 = @(x) -4253*x.^3+2262*x.^2-400.2*x+23.31;

    flashfit = @(x) p1(x).*(x<=0.046699) + ...
                p2(x).*(x> 0.046699).*(x<=0.1185) +...
                p3(x).*(x>0.1185);
            
    % Compute l2 error of this trial with polynomial
    %  Find times in coordinates of model that are in domain of fit
    flag_t = logical((taxis >= (t_interval(1)-tshift)).*...
                 (taxis <= (t_interval(2))-tshift));

    % Note experimental trace placed a flash at 5 ms
    tpts = taxis(flag_t);
    %  To evaluate fit go back to domain of fit
    tpts = tpts + tshift;

    simpts = J_tot(1)-J_tot; 
    simpts = (1e12)*simpts;
    simpts = simpts(flag_t);

    flashpts = flashfit(tpts);

    % Store
    functionals(7) = sqrt(...
                         sum(...
                             (flashpts - simpts).^2 ...
                            )...
                     );
    case false
        functionals(7) = NaN;
end

if functionals(7) > 1e3
    functionals(7) = 1e3;
end

%% Tpeak for J_tot
[~,tpeak] = max(J_drop);
tpeak = taxis(tpeak);

if tpeak > 0.25
    tpeak = 0.25;
end

% Store 
functionals(8) = tpeak;

%% Iwidth
% Create a flag to say whether overshoot is real
%   flag_over is used at other points in script
Ipts = J_tot(1) - J_tot;
Imin = min(Ipts);
flag_over = ( Imin < 0 )&&( Imin < Ipts(end) ); 

%   Find last time point before drop was first negative
Ipts = J_tot(1) - J_tot;
pos = find(Ipts < 0,1);
if (flag_over)&&(~isempty(pos))
    delta = taxis(pos-1);
else
    delta = taxis(max(size(taxis)));
end

% Store
functionals(9) = delta;

%% Jdark
jd = J_tot(1)*1e12; % Convert A to pA

if jd < 0
    jd = 0;
elseif jd > 48
    jd = 48;
end

% Store
functionals(10) = jd;

%% Overshoot of Jdrop
jover = abs(min(min(J_drop),0));

if jover > 6.5
    jover = 6.5;
end

% Store
functionals(11) = jover*flag_over;

%% Indicator Ipeak
thresh = 13.3; % pA (one-half saturating value in exp)
peak = functionals(2);

functionals(12) = (peak >= thresh);

%% Indicator L2 Fit
% Check that jdark is within 10% and the average DROP error is within 
%   thresh compared to Ingram et al
thresh = .5;
jd = functionals(10);

% Do the calculation assuming simpts etc are defined
try
    simdrop   = simpts/functionals(10);
    flashdrop = flashpts/26.6;
    
    % Average error across time points
    %   Since the drop is mostly between 0 and 1 (and in Ingram et al is
    %   always between 0 and 1) we're saying the difference is more or less
    %   bounded by 0 and 1, hence thresh is a number in that interval.
    errdrop = sqrt(...
                   sum(...
                       (flashdrop-simdrop).^2 ...
                      )...
                      /length(simdrop)...
                  );
    
    functionals(13) = (errdrop <= thresh)&(jd >= .9*26.6)&(jd <= 1.1*26.6);
    
catch
    functionals(13) = NaN;
end
%% Indicator of tpeak
thresh = .01; % +- 10 ms rel Ingram et al for when time of peak drop occurs
tpeak = functionals(8);

functionals(14) = (tpeak >= .041 - thresh)&(tpeak <= .041 + thresh);

%% Indicator of Jdark
thresh = .5; % +- factor of magnitude of jdark rel Ingram et al
jd = functionals(10); % Convert A to pA

functionals(15) = (jd >= thresh*26.6)&(jd <= (1+thresh)*26.5);

%% Indicator of overshoot
jover = functionals(11);

functionals(16) = (jover > 0);
end

