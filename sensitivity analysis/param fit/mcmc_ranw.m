function [smp,prange,threshold,p,gamma]...
              = mcmc_ranw(flag_params,smp,pointer,fitscore)
%Propose a new sample based at current in random walk fashion
%   flag_params is 1 x n_datamat with true indicating which entries of
%    sample should be drawn from the below probability distributions
%   smp is formatted as datamat and it says where to base point for the
%    Metropolis Hastings transition kernel
%   pointer is exported from datamat
%   smp as output is the input but updated to reflect the drawn samples
%   prange gives the *desired* range of parameters in two column format. It
%    is used to define the standard deviations in the randw's normal
%    distributions
%   threshold defines value so that if fitscore is above we have a chance
%    propose from a uniform rather than Gaussian ranw

%% Proposal kernel parameters
% Fitscore threshold above which we may draw a uniform sample
threshold = Inf;

% Probability of drawing uniform sample if above threshold
%  1/p is mean of geometric distribution
p = 1/5; % On average will draw uniform sample by 20th try

% Ranw's std deviation as percentage of desired interval length
gamma = .01; % 2 std dev's is 14% of given prange 

%% Initialize parameters
% So that mc_sample treats j_cG_max and j_ex_sat in standard way
%  Normalize out the 1/(B_Ca*F) factor included by datamat in case one of
%  these parameters is not varied so then final lines of code don't double
%  scale by this factor
B = pointer(12);
smp(B+1) = smp(B+1)*...
                   (smp(pointer(5)+2)*smp(pointer(4)+2)); % j_cg_max
                 
B = pointer(13);
smp(B+1) = smp(B+1)*...
                   (smp(pointer(5)+2)*smp(pointer(4)+2)); % j_ex_sat

% Count how many values are varied
n_smp = nnz(flag_params);

%% Continue if doing more than flag_params adjust
if nargin == 1
    smp = NaN;
    prange = NaN;
else
%% Give preferred parameter ranges
prange = zeros(n_smp,2);

% Geometric params
prange(1:12,:) = ...
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
prange(B+1:B+4,:) = [25 25; ... % n_sez
                   .6/1.5 .6/1.5; ... % taglia
                   1e-5 1e-5; ... % tol_R
                   1e-3 1e-3; ... % tol_angle
                  ];
              
% Time discretization params
B = pointer(2);
prange(B+1:B+10,:) = [0 0; ... % method_cyto(1)
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
prange(B+1,:) = [1e-8 1e-8]; % tol_stat

% Physical constants
B = pointer(4);
prange(B+1:B+2,:) = [6.022e23 6.022e23; ... % N_Av
                     96500/1e21 96500/1e21; ... % F
                    ];
                
% Cytoplasmic buffer
B = pointer(5);
prange(B+1:B+2,:) = [...
                     1 2; ... % B_cG
                     10 30; ... % B_Ca
                    ];
                
% R_st biochemistry
B = pointer(6);
prange(B+1:B+4,:) = [...
                     30 330;... % nu_RG
                     1 2; ... % D_R_st
                     1 200; ... % k_R
                     99.5e3 100.5e3; ... R_sigma
                    ];
                
% G_st biochemistry
B = pointer(7);
prange(B+1:B+3,:) = [...
                     .05 1; ... % k_GE
                     1.2 3.2; ... % D_G_st
                     500 1500; ... % G_sigma
                    ];
                
% PDE_st biochemistry (free)
B = pointer(8);
prange(B+1:B+6,:) = [...
                     .8 1.6; ... % D_E_st
                     5 150; ... % k_E (Majumder15 -> 2.3xShen )
                     10 100; ... % PDE_sigma
                     1 1e3; ... % Beta_dark
                     10e-6 10e-6; ... % K_m
                     8.3e-6 8.3e-6; ... % k_st
                    ];
                
% PDE_st biochemistry (derived)
B = pointer(9);
prange(B+1:B+5,:) = [...
                     1 1; ... % E_sigma
                     1 1; ... % E_vol
                     1 1; ... % k_hyd
                     .16 1.5; ... % kcat_DIV_Km
                     1 1; ... % k_st
                    ];
                
% cG biochemistry
B = pointer(10);
prange(B+1:B+2,:) = [...
                     2 4; ... % u_tent
                     50 196; ... % kk_u
                    ];
                
% Ca biochemistry
B = pointer(11);
prange(B+1:B+2,:) = [...
                     .2 .67; ... % v_tent
                     15 15; ... % kk_v
                    ];
                
% cG-gated current
B = pointer(12);
prange(B+1:B+4,:) = [...
                     1000e-12 5000e-12; ... % j_cg_max
                     2.5 3.5; ... % m_cG
                     10 30; ... % K_cG
                     .2 .35; ... % f_Ca
                    ];
                
% Exchanger current
B = pointer(13);
prange(B+1:B+2,:) = [...
                     1e-12 10e-12; ... % j_ex_sat
                     .02 5 ... % K_ex
                    ];
                
% Cyclase biochemistry
B = pointer(14);
prange(B+1:B+4,:) = [...
                     50 500; ... % alpha_max
                     2 20; ... % alpha_max/alpha_min
                     2 2.5; ... % m_cyc
                     130e-3 140e-3; ... % K_cyc
                    ];
                
% Illumination 
B = pointer(15);
prange(B+1:B+2,:) = [...
                     355 355; ... % n_Rst0
                     0 0; ... % rate
                    ];
% Restart param space-holders
B = pointer(16);
prange(B+1:B+2,:) = [...
                     0 0; ... % flag_restart
                     0 0; ... % rst_index
                    ];
                
%% Implement proposal function
if (fitscore >= threshold)&&(rand(1) <= p)
    % Draw uniform samples from prange
    unifs = prange(flag_params,1) + ...
            (prange(flag_params,2) - prange(flag_params,1)).*rand(n_smp,1);
    smp(flag_params) = unifs;
    
else
    % RanW by normal
    % Revert alpha_min position to the alpha_max/alpha_min value so ranw
    %  as intended
    B = pointer(14);
    smp(B+2) = smp(B+1)/smp(B+2);
    
    %  Standard deviation is gamma fraction of interval length
    stddev = gamma*(prange(:,2) - prange(:,1));

    % Sample normal distribution with desired standard deviations
    normals = stddev(flag_params).*randn(n_smp,1);

    % Produce next sample
    smp(flag_params) = smp(flag_params) + normals;
end

%% Set dep params by ind param values
smp(4) = fcosgamma0(smp(1),smp(2),smp(3)); % cosgamma0

B = pointer(9);
smp(B+1) = fE_sigma(smp(pointer(8) + 3)); % E_sigma
smp(B+2) = fE_vol(smp(pointer(8) + 3),smp(8),smp(7)); % E_vol
smp(B+3) = fk_hyd(smp(8),smp(7),smp(pointer(8)+4),smp(pointer(8)+3)); % k_hyd
smp(B+5) = fk_st(smp(B+4),smp(pointer(5)+1)); % k_st

B = pointer(12);
smp(B+1) = smp(B+1)/...
                   (smp(pointer(5)+2)*smp(pointer(4)+2)); % j_cg_max
                 
B = pointer(13); % Last block completed
smp(B+1) = smp(B+1)/...
                   (smp(pointer(5)+2)*smp(pointer(4)+2)); % j_ex_sat

B = pointer(14); % Last block completed
smp(B+2) = smp(B+1)/( smp(B+2) ); % alpha_min

end

end

