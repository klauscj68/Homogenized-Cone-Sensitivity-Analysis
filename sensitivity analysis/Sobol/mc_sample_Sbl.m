function [smp] = mc_sample_Sbl(flag_params,smp,pointer)
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

% So that we have true independent parameters for Sobol analysis,
%  check normalization of amin to [0,1]
B = pointer(14);
if (smp(B+2)<0)||(smp(B+2)>1)
    warning('The passed amin value is inadmissible'); 
end

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

%% Geometric param r.v.'s
Distr{1} = @(x) (1-x)*.585+x*.615; % R_b
Distr{2} = @(x) (1-x)*.38+x*.42; % R_t
Distr{3} = @(x) (1-x)*12.7+x*14.1; % H
Distr{4} = @(x) 1; % cosgamma0
Distr{5} = @(x) 0; % theta_in
Distr{6} = @(x) (1-x)*.8*pi+x*1.2*pi; % theta_fin
Distr{7} = @(x) (1-x)*.0164+x*.0172; % epsilon_0
Distr{8} = @(x) (1-x)*.61+x*(.71); % nu
Distr{9} = @(x) (1-x)*.8+x*1.2; % sigma
Distr{10} = @(x) 1; %flag_ch
Distr{11} = @(x) 0; %flag_ch
Distr{12} = @(x) 0; %flag_ch

%% Mesh discretization r.v.'s
B = pointer(1); % Last block completed
Distr{B+1} = @(x) 25; % n_sez
Distr{B+2} = @(x) 0.4; % taglia
Distr{B+3} = @(x) 2.4e-4; % tol_R
Distr{B+4} = @(x) 6.2832e-3; % tol_angle

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
Distr{B+9} = @(x) 1350; % n_step_t
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
Distr{B+1} = @(x) (1-x)*1+x*2; % B_cG
Distr{B+2} = @(x) (1-x)*10+x*30; % B_Ca

%% R_st biochemistry r.v.'s
B = pointer(6); % Last block completed
Distr{B+1} = @(x) (1-x)*30+x*330; % nu_RG
Distr{B+2} = @(x) (1-x)*1+2*x; %D_R_st
Distr{B+3} = @(x) (1-x)*1+x*200; % k_R 
Distr{B+4} = @(x) (1-x)*8e3+x*(1.2e4); % R_sigma

%% G_st biochemistry r.v.'s
B = pointer(7); % Last block completed
Distr{B+1} = @(x) (1-x)*.05+1*x; % k_GE
Distr{B+2} = @(x) (1-x)*1.1+x*3.2; % D_G_st
Distr{B+3} = @(x) (1-x)*200+x*1500; % G_sigma

%% PDE_st biochemistry r.v.'s (free params)
B = pointer(8); % Last block completed
Distr{B+1} = @(x) (1-x)*.8+x*1.6; % D_E_st
Distr{B+2} = @(x) (1-x)*5+x*150; % k_E
Distr{B+3} = @(x) (1-x)*10+x*120; % PDE_sigma
Distr{B+4} = @(x) (1-x)*1+x*1e3; % Beta_dark
Distr{B+5} = @(x) (.95+.1*x)*10e-6; % K_m
Distr{B+6} = @(x) (.95+.1*x)*4.1756e-5; % k_cat

%% PDE_st biochemistry r.v.'s (derived params)
% Except for kcat_DIV_Km these don't matter because they get overwritten in
% final section
B = pointer(9); % Last block completed
Distr{B+1} = @(x) 1; % E_sigma
Distr{B+2} = @(x) 1; % E_vol
Distr{B+3} = @(x) 1; % k_hyd
Distr{B+4} = @(x) ((1-x)*95+x*905)/(6.022e2); %kcat_DIV_Km for E* half-unit
Distr{B+5} = @(x) 1; % k_st

%% cG biochemistry r.v.'s
B = pointer(10); % Last block completed
Distr{B+1} = @(x) (.95+.1*x)*2; % u_tent
Distr{B+2} = @(x) (1-x)*50+x*196; % kk_u

%% Ca biochemistry r.v.'s
B = pointer(11); % Last block completed
Distr{B+1} = @(x) (.95+.1*x)*.5; % v_tent
Distr{B+2} = @(x) (1-x)*12+x*18; % kk_v

%% cG-gated current r.v.'s
B = pointer(12); % Last block completed
Distr{B+1} = @(x) ((1-x)*1e3+x*5e3)*1e-12; % j_cg_max
Distr{B+2} = @(x) (1-x)*2+x*3.5; % m_cG
Distr{B+3} = @(x) (1-x)*10+x*30; % K_cG
Distr{B+4} = @(x) (1-x)*.2+x*.35; % f_Ca

%% Exchanger current r.v.'s
B = pointer(13); % Last block completed
Distr{B+1} = @(x) ((1-x)*1+x*10)*1e-12; % j_ex_sat
Distr{B+2} = @(x) (1-x)*.02+x*5; % K_ex

%% Cyclase biochemistry r.v.'s
B = pointer(14); % Last block completed
Distr{B+1} = @(x) (1-x)*50+x*500; % alpha_max
Distr{B+2} = @(x) x; % COV'd alpha_min
Distr{B+3} = @(x) (1-x)*2+x*2.5; % m_cyc
Distr{B+4} = @(x) (1-x)*.130+x*.140; %K_cyc

%% Illumination r.v.'s
B = pointer(15); % Last block completed
Distr{B+1} = @(x) 940; % n_Rst0
Distr{B+2} = @(x) 0; %rate

%% Restart param space-holders
B = pointer(16); % Last block completed
Distr{B+1} = @(x) 0; % flag_restart
Distr{B+2} = @(x) 0; % rst_index

%% Map uniform samples by r.v.'s like the maps above
update = find(flag_params);
n_update = length(update);

unif = rand(n_update,1);

for i=1:n_update
    pos = update(i);
    smp(pos) = Distr{pos}(unif(i));
end

%% Set dep params by ind param values
smp(4) = fcosgamma0(smp(1),smp(2),smp(3)); % cosgamma0

B = pointer(9);
smp(B+1) = fE_sigma(smp(pointer(8) + 3)); % E_sigma
smp(B+2) = fE_vol(smp(pointer(8) + 3),smp(8),smp(7)); % E_vol
smp(B+3) = fk_hyd(smp(8),smp(7),smp(pointer(8)+4),smp(pointer(8)+3)); % k_hyd
smp(B+5) = fk_st(smp(B+4),smp(pointer(5)+1)); % k_st
      
% Return j_cg_max and j_ex_sat to their normalized form
B = pointer(12);
smp(B+1) = smp(B+1)/...
                   (smp(pointer(5)+2)*smp(pointer(4)+2)); % j_cg_max
                 
B = pointer(13); % Last block completed
smp(B+1) = smp(B+1)/...
                   (smp(pointer(5)+2)*smp(pointer(4)+2)); % j_ex_sat

end

