function [flag_accept,flag_nstay] = mcmc_alpha(fitscore1,fitscore2,...
                                    base,propose,pointer,...
                                    prange,flag_params,...
                                    threshold,p,pr_gamma)
%Metropolis-Hastings acceptance random variable for symmetric transition 
%   fitscore1, fitscore2 are fitness values for the base proposed sample
%    respectively.
%   propose is the datamat for the proposed sample; pointer from data_set
%   base is the same for the base point parameter set
%   prange is the desired range of parameters we'd prefer to sample, seen
%    through an exponential decay outside of the area. 
%   flag_params is true or false to indicate which parameters are being
%    varied in the walk. Needed so if doing a uniform draw, we can find
%    normalizing factor for density computation
%   threshold is mcmc_ranw output which defines the fitness level above
%    which there's a chance we do a uniform draw rather than a ran walk
%   p defines the probability of transitioning to uniform draw from prange
%    if the fitness score is above threshold
%   pr_gamma is the percentage of prange interval used for standard
%    deviation in proposal function's ranw
%   flag_accept is false or true based on whether the candidate has been
%   randomly accepted
%   flag_nstay is a dummy var, simply negating flag_accept
%   Code works by proposing a desired stationary distribution density pi
%    ~ \Prod exp(-gamma*dist(parami,prange(:,i))/length(prange(:,i))*...
%                (dist(...) >= 0)*...
%            (1/fitness)^beta
%    beta is chosen so that if proposed fitness score is 1.25 the base
%     value it is only accepted 5% of the time when in the desired interval
%    gamma is chosen so that if 200% the length away from desired interval
%     you are only 25% as likely to accept as before 

%% Initialize values
beta = 13.425; % log(.05)/log(1/1.25), ie (1.25)^beta = .05
gamma = .6931; % exp(-2*gamma) = .25

n_smp = max(size(base));

% Convert alpha_min entry into alpha_max/alpha_min value
B = pointer(14);
propose(B+2) = propose(B+1)/propose(B+2);

%% Set mandatory constraints on parameters if any
% Setting flag_cnst to true means the parameter must fall in range below
%  NOTE: wellposed_ss also sets absolute constraints on the parameters in
%        that the cG_d, Ca_d, j_dark mandated ranges are given. Also,
%        params are required to be positive there
flag_cnst = false(n_smp,1);

B = pointer(5);
flag_cnst(B+1) = true; % B_cG
flag_cnst(B+2) = true; % B_Ca

B = pointer(6);
flag_cnst(B+1) = true; % nu_RG
flag_cnst(B+3) = true; % k_R

B = pointer(7);
flag_cnst(B+1) = true; % k_GE
flag_cnst(B+3) = true; % G_sigma

B = pointer(8);
flag_cnst(B+2) = true; % k_E
flag_cnst(B+3) = true; % PDE_sigma
flag_cnst(B+4) = true; % Beta_dark

B = pointer(9);
flag_cnst(B+4) = true; % kcat/Km

B = pointer(12);
flag_cnst(B+1) = true; % j_cG_max
flag_cnst(B+2) = true; % m_cG
flag_cnst(B+3) = true; % K_cG
flag_cnst(B+4) = true; % f_Ca

B = pointer(13);
flag_cnst(B+1) = true; % j_ex_sat
flag_cnst(B+2) = true; % K_ex

B = pointer(14);
flag_cnst(B+1) = true; % amax
flag_cnst(B+2) = true; % amax/amin
flag_cnst(B+3) = true; % m_cyc
flag_cnst(B+4) = true; % K_cyc

% Give preferred parameter ranges
Crange = zeros(n_smp,2);

% Geometric params
Crange(1:12,:) = ...
                 [3 3; ... % R_b
                  1 1; ... % R_t
                  15 15; ... % H
                  1 1; ... % cosgamma0
                  0 0; ... % theta_in
                  pi pi; ... % theta_fin
                  .015 .015; ... % epsilon_0
                  1 1; ... % nu
                  1 1; ... % sigma
                  1 1; ... % flag_ch
                  0 0; ... % flag_ch
                  0 0; ... % flag_ch
                 ];
             
% Mesh params
B = pointer(1); 
Crange(B+1:B+4,:) = [25 25; ... % n_sez
                   .6/1.5 .6/1.5; ... % taglia
                   1e-5 1e-5; ... % tol_R
                   1e-3 1e-3; ... % tol_angle
                  ];
              
% Time discretization params
B = pointer(2);
Crange(B+1:B+10,:) = [0 0; ... % method_cyto(1)
                      25 25; ... % method_cyto(2)
                      5 5; ... % method_cyto(3)
                      .5 .5; ... % theta
                      1 1; ... % alpha
                      1e-3 1e-3; ... % tol_fix
                      0 0; ... % norma_inf
                      .5 .5; ... % t_fin
                      3*450 3*450; ... % n_step_t
                      1 1; ... % downsample
                     ];
                 
% Steady state solver
B = pointer(3);
Crange(B+1,:) = [1e-8 1e-8]; % tol_stat

% Physical constants
B = pointer(4);
Crange(B+1:B+2,:) = [6.022e23 6.022e23; ... % N_Av
                     96500/1e21 96500/1e21; ... % F
                    ];
                
% Cytoplasmic buffer
B = pointer(5);
Crange(B+1:B+2,:) = [...
                     1 2; ... % B_cG
                     10 30; ... % B_Ca
                    ];
                
% R_st biochemistry
B = pointer(6);
Crange(B+1:B+4,:) = [...
                     30 300;... % nu_RG
                     1 2; ... % D_R_st
                     1 200; ... % k_R
                     10^5 10^5; ... R_sigma
                    ];
                
% G_st biochemistry
B = pointer(7);
Crange(B+1:B+3,:) = [...
                     .05 1; ... % k_GE
                     1.2 3.2; ... % D_G_st
                     200 1500; ... % G_sigma
                    ];
                
% PDE_st biochemistry (free)
B = pointer(8);
Crange(B+1:B+6,:) = [...
                     .8 1.6; ... % D_E_st
                     5 150; ... % k_E (Majumder15+Shen -> 1-7x rod)
                     10 120; ... % PDE_sigma
                     1 1000; ... % Beta_dark
                     10e-6 10e-6; ... % K_m
                     8.3e-6 8.3e-6; ... % k_st
                    ];
                
% PDE_st biochemistry (derived)
B = pointer(9);
Crange(B+1:B+5,:) = [...
                     1 1; ... % E_sigma
                     1 1; ... % E_vol
                     1 1; ... % k_hyd
                     0.1577 1.5028;...%.5*830e-3 5*830e-3; ... % kcat_DIV_Km
                     1 1; ... % k_st
                    ];
                
% cG biochemistry
B = pointer(10);
Crange(B+1:B+2,:) = [...
                     2 4; ... % u_tent
                     50 196; ... % kk_u
                    ];
                
% Ca biochemistry
B = pointer(11);
Crange(B+1:B+2,:) = [...
                     .2 .67; ... % v_tent
                     15 15; ... % kk_v
                    ];
                
% cG-gated current
B = pointer(12);
Crange(B+1:B+4,:) = [...
                     1000e-12 5000e-12; ... % j_cg_max
                     2 3.5; ... % m_cG
                     10 30; ...%26 504; ... % K_cG (Haynes Yau 1990 44uM)
                     .2 .35;...%.28 .4; ... % f_Ca (Ohyama .34+-.06 range)
                    ];
                
% Exchanger current
B = pointer(13);
Crange(B+1:B+2,:) = [...
                     1e-12 10e-12; ... % j_ex_sat
                     .02 5 ... % K_ex
                    ];
                
% Cyclase biochemistry
B = pointer(14);
Crange(B+1:B+4,:) = [...
                     50 500; ... % alpha_max
                     2 20; ... % alpha_max/alpha_min (Peshenko et al, 28x but Kawamura less)
                     2 2.5; ... % m_cyc
                     130e-3 140e-3; ... % K_cyc
                    ];

% Illumination 
B = pointer(15);
Crange(B+1:B+2,:) = [...
                     372 372; ... % n_Rst0
                     0 0; ... % rate
                    ];
%% Renormalize j's                
% So that mc_sample treats j_cG_max and j_ex_sat in standard way
%  Normalize out the 1/(B_Ca*F) factor included by datamat in case one of
%  these parameters is not varied so then final lines of code don't double
%  scale by this factor
B = pointer(12);
propose(B+1) = propose(B+1)*...
                   (propose(pointer(5)+2)*propose(pointer(4)+2)); % j_cg_max
               
base(B+1) = base(B+1)*...
                   (base(pointer(5)+2)*base(pointer(4)+2)); % j_cg_max
                 
B = pointer(13);
propose(B+1) = propose(B+1)*...
                   (propose(pointer(5)+2)*propose(pointer(4)+2)); % j_ex_sat
base(B+1) = base(B+1)*...
                   (base(pointer(5)+2)*base(pointer(4)+2));
                
                
%% Compute acceptance factor
% The constraint ranges part of MH-ratio
if nnz(flag_cnst)~=0
    % Check if anywhere a value fell outside of range
    ram = (propose(flag_cnst) < Crange(flag_cnst,1))|...
           (propose(flag_cnst) > Crange(flag_cnst,2));
       
       if (real(fitscore2) >= 250)||(nnz(ram) ~= 0)
           flag_accept = false;
           flag_nstay = ~flag_accept;
           return
       end
end

% The fitness part of the Metropolis-Hastings ratio
alpha = (fitscore1/fitscore2)^beta;

% The exponential decay support part
idx = find(flag_params);
decay1 = 1;
decay2 = 1;
for i=1:max(size(idx))
    % Base point factor
    pdist = reldist(base(idx(i)),prange(idx(i),:));
    decay1 = exp(-gamma*pdist)*decay1;
    
    % Proposed point factor
    pdist = reldist(propose(idx(i)),prange(idx(i),:));
    decay2 = exp(-gamma*pdist)*decay2;
end
decay = decay2/decay1;

% The asymmetric uniform-normal proposal part
if ((fitscore1 <= threshold)&&(fitscore2 <= threshold))||...
        ((fitscore1 > threshold)&&(fitscore2 > threshold))
    % Points below are at same threshold; their proposal densities are same
    asym = 1;
elseif fitscore1 <= threshold
    asym = 1;
    for i=1:max(size(idx))
        % Compute normalizing length of prange interval and associated std
        vol = prange(idx(i),2)-prange(idx(i),1);
        sigstd = pr_gamma*vol;
        asym = asym*(...
                     1-p + p/vol*(...
                            sigstd*sqrt(2*pi)*...
                            exp(.5*((base(idx(i)) - propose(idx(i)))/sigstd)^2)...
                                 )*...
                            (base(idx(i))>=prange(idx(i),1))*(base(idx(i))<=prange(idx(i),2))...
                     );
    end
else
    asym = 1;
    for i=1:max(size(idx))
        % Compute normalizing length of prange interval and associated std
        vol = prange(idx(i),2)-prange(idx(i),1);
        sigstd = pr_gamma*vol;
        asym = asym*(...
                     1-p + p/vol*(...
                            sigstd*sqrt(2*pi)*...
                            exp(.5*((base(idx(i)) - propose(idx(i)))/sigstd)^2)...
                                 )*...
                            (propose(idx(i))>=prange(idx(i),1))*(propose(idx(i))<=prange(idx(i),2))...
                     )^(-1);
    end
end

% Metropolis Hastings choice of alpha
alpha = min(decay*alpha*asym,1);

%% Decide if sample accepted
if rand(1) <= alpha
    flag_accept = true;
else
    flag_accept = false;
end
flag_nstay = ~flag_accept;
end

%% Relative distance function
function [rdist] = reldist(val,pint)
% Compute relative distance of val from interval
% Make sure pint ordered least to greatest
if pint(2) <= pint(1)
    pint = [pint(2) pint(1)];
end

length = pint(2) - pint(1);

if val < pint(1)
    rdist = (pint(1)-val)/length;
elseif val > pint(2)
    rdist = (val - pint(2))/length;
else
    rdist = 0;
end

end
