function [smp] = mc_sample(flag_params,smp,pointer)
%Sample desired parameters like prescribed random variables
%   flag_params is 1 x n_datamat with true indicating which entries of
%    sample should be drawn from the below probability distributions
%   smp is essentially values like datamat with role to say what values
%    that arent to be sampled should take
%   pointer is exported from datamat
%   smp as output is the input but updated to reflect the drawn samples

%% Initialize parameters
n_smp = max(size(smp));
Distr = cell(n_smp,1);

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



% Don't vary the parameters which aren't free/are determined
flag_params(4) = false; % cosgamma0
flag_params(5) = false; % theta_in

B = pointer(1);
flag_params(B+3:B+4) = false; % tol_angle, tol_R

B = pointer(2);
flag_params(B+1:B+7) = false; % method_cyto, theta, alpha, tol_fix, norma
flag_params(B+10) = false; % downsample

B = pointer(3);
flag_params(B+1) = false; % tol_stat

B = pointer(4);
flag_params(B+1:B+2) = false; % N_Av,F

B = pointer(8);
flag_params(B+5:B+6) = false; % K_m, k_cat

B = pointer(9);
flag_params(B+1:B+3) = false; % E_sigma, E_vol, k_hyd
flag_params(B+5) = false; % k_st

B = pointer(15);
flag_params(B+1:B+2) = false; % n_Rst0, rate

B = pointer(16);
flag_params(B+1:B+2) = false; % flag_restart, rst_index

%% Take uniform [0,1] r.v. samples to later map to desired 
unif = rand(n_smp,1);

%% Geometric param r.v.'s
Distr{1} = @(x) .6; % R_b
Distr{2} = @(x) .4; % R_t
Distr{3} = @(x) 13.4; % H
Distr{4} = @(x) 1; % cosgamma0
Distr{5} = @(x) 0; % theta_in
Distr{6} = @(x) pi; % theta_fin
Distr{7} = @(x) .0134 + x*(.0168-.0134); % epsilon_0
Distr{8} = @(x) 11/16.8+x*(1-11/16.8); % nu
Distr{9} = @(x) 1; % sigma
Distr{10} = @(x) 1; %flag_ch
Distr{11} = @(x) 0; %flag_ch
Distr{12} = @(x) 0; %flag_ch

%% Mesh discretization r.v.'s
B = pointer(1); % Last block completed
Distr{B+1} = @(x) 50; % n_sez
Distr{B+2} = @(x) .6/3.5; % taglia
Distr{B+3} = @(x) 1e-5; % tol_R
Distr{B+4} = @(x) 1e-3; % tol_angle

%% Time discretization r.v.'s
B = pointer(2); % Last block completed
Distr{B+1} = @(x) 0; % method_cyto(1)
Distr{B+2} = @(x) 25; % method_cyto(2)
Distr{B+3} = @(x) 5; % method_cyto(3)
Distr{B+4} = @(x) .5; % theta
Distr{B+5} = @(x) 1; % alpha
Distr{B+6} = @(x) 1e-3; % tol_fix
Distr{B+7} = @(x) 0; % norma_inf
Distr{B+8} = @(x) .5; % t_fin
Distr{B+9} = @(x) 3*450; % n_step_t
Distr{B+10} = @(x) 1; % downsample
 
%% Steady state solver r.v.'s
B = pointer(3); % Last block completed
Distr{B+1} = @(x) 1e-8; % tol_stat

%% Physical constants r.v.'s
B = pointer(4); % Last block completed
Distr{B+1} = @(x) 6.022e23; % N_Av
Distr{B+2} = @(x) 96500/1e21; % F

%% Cytoplasmic buffer r.v.'s
B = pointer(5); % Last block completed
Distr{B+1} = @(x) 1 + x*(2 - 1); % B_cG
Distr{B+2} = @(x) 10 + x*(44 - 10); % B_Ca

%% R_st biochemistry r.v.'s
B = pointer(6); % Last block completed
Distr{B+1} = @(x) 30+x*(330-30); % nu_RG
Distr{B+2} = @(x) 1 + x*(2-1); %D_R_st
Distr{B+3} = @(x) 25/4+x*(65-25/4); % k_R 
Distr{B+4} = @(x) 10^5; % R_sigma

%% G_st biochemistry r.v.'s
B = pointer(7); % Last block completed
Distr{B+1} = @(x) 1; % k_GE
Distr{B+2} = @(x) 1.2 + x*(3.2-1.2); % D_G_st
Distr{B+3} = @(x) 2000 + x*(3000-2000); % G_sigma

%% PDE_st biochemistry r.v.'s (free params)
B = pointer(8); % Last block completed
Distr{B+1} = @(x) .8 + x*(1.6-.8); % D_E_st
Distr{B+2} = @(x) 18.5/4 + x*(55-18.5/4); % k_E
Distr{B+3} = @(x) 500 + x*(1000 - 500); % PDE_sigma
Distr{B+4} = @(x) 10 + x*(100-10); % Beta_dark
Distr{B+5} = @(x) 10e-6; % K_m
Distr{B+6} = @(x) 8.3e-6; % k_cat

%% PDE_st biochemistry r.v.'s (derived params)
% Except for kcat_DIV_Km these don't matter because they get overwritten in
% final section
B = pointer(9); % Last block completed
Distr{B+1} = @(x) 1; % E_sigma
Distr{B+2} = @(x) 1; % E_vol
Distr{B+3} = @(x) 1; % k_hyd
Distr{B+4} = @(x) .5*830e-3+x*(5*830e-3-.5*830e-3); %kcat_DIV_Km
Distr{B+5} = @(x) 1; % k_st

%% cG biochemistry r.v.'s
B = pointer(10); % Last block completed
Distr{B+1} = @(x) 2+x*(4-2); % u_tent
Distr{B+2} = @(x) 50 + x*(196-50); % kk_u

%% Ca biochemistry r.v.'s
B = pointer(11); % Last block completed
Distr{B+1} = @(x) .2 + x*(.67-.2); % v_tent
Distr{B+2} = @(x) 15; % kk_v

%% cG-gated current r.v.'s
B = pointer(12); % Last block completed
Distr{B+1} = @(x) (2000e-12+x*(8000e-12 - 2000e-12)); % j_cg_max
Distr{B+2} = @(x) 2.5 + x*(3-2.5); % m_cG
Distr{B+3} = @(x) 10+x*(30-10); % K_cG
Distr{B+4} = @(x) .2+x*(.35-.2); % f_Ca

%% Exchanger current r.v.'s
B = pointer(13); % Last block completed
Distr{B+1} = @(x) (1e-12 + x*(10e-12 - 1e-12)); % j_ex_sat
Distr{B+2} = @(x) .9 + x*(1.6-.9); % K_ex

%% Cyclase biochemistry r.v.'s
B = pointer(14); % Last block completed
Distr{B+1} = @(x) 76.5 + x*(800-76.5); % alpha_max
Distr{B+2} = @(x) ( (6.7+x*(13.9-6.7)) ); % alpha_max/alpha_min
Distr{B+3} = @(x) (1.6 + x*(3-1.6)); % m_cyc
Distr{B+4} = @(x) 73e-3 + x*(400e-3 - 73e-3); %K_cyc

%% Illumination r.v.'s
B = pointer(15); % Last block completed
Distr{B+1} = @(x) 940; % n_Rst0
Distr{B+2} = @(x) 0; %rate

%% Restart param space-holders
B = pointer(16); % Last block completed
Distr{B+1} = @(x) 0; % flag_restart
Distr{B+2} = @(x) 0; % rst_index

%% Map uniform samples by r.v.'s above
update = find(flag_params);
n_update = max(size(update));

for i=1:n_update
    pos = update(i);
    smp(pos) = Distr{pos}(unif(pos));
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

