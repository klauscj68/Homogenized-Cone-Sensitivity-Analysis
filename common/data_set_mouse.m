function [datamat,pointer,...
          fnames] = data_set()
%Where parameters are specified for the bright light cone simulation.
%   datac is a n x 1 array whose entries group the different parameters by
%    type. We use array because this in tandem with data.m make it easier
%    to write parameter sweeps using Matlab loops.
%   pointer is a #sections by 1 array whose entry value at i is the last
%    index of datac belonging to i^th section. 
%   fnames is the character names of the .mat files needed to restart the
%   simulation

%% Initialize datac and pointer
datamat =   [];
pointer = [];

%% Geometric Parameters
% Radius at the bottom
R_b = .6; % um

% Radius at the top
R_t = .4; % um

% Cone Height
H = 13.4; % um

% Aperture angle of cone
h = abs(R_t*H/(R_b-R_t));
if isinf(h)
    cosgamma0 = 1;
else
    cosgamma0 = h/sqrt(h^2+R_t^2);
end

% Connecting sliver angular interval
theta_in  = 0; 
theta_fin = pi;

% Disc thickness
%  Formula is H/((n_ch-1)*2) where n_ch is the number of chambers
epsilon_0 = .0168; % um

% Ratio between interdiscal and disc thickness
nu = 11/16.8;

% Ratio between outer-shell and disc thickness
sigma = 15/15;

% Location of channels
%  In order, says if channels reside at sliver, discs, lateral folds
flag_ch = [true false false]; 

% Record into datac
datamat = [datamat;...
         R_b;...
         R_t;...
         H;...
         cosgamma0;...
         theta_in;...
         theta_fin;...
         epsilon_0;...
         nu;...
         sigma;...
         flag_ch(:)];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Space Discretization
% Number of z-horizontal cross sections
n_sez = 25;

% Target maximum edge length of triangles in cross section
taglia = R_b/1.5 ;

% Linear tolerance of grid for searching nodes
%  Nodes are considered same if their distance is less than tol_R
tol_R = taglia/1000; %\um

% Angular tolerance of grid for searching nodes
tol_angle = 2*pi/1000;

% Record into datac
datamat = [datamat;...
         n_sez;...
         taglia;...
         tol_R;...
         tol_angle];

% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Time Integration
% Integration method
% method_cyto(1) indicates the type of solver
% 0 means Matlab built in solver ODE45.
%  If 0, then nextentries of method_cyto(2:3) are not used
% 1 means solve using theta method
% method_cyto(2) says how many consecutive samples are used to assess
%  Cauchy convergence of the fixed point solution for S_theta
% method_cyto(3) says how many samples to average together before assessing
%  Cauchy convergence
method_cyto = [0 25 5];
assert( floor(method_cyto(2)/method_cyto(3)) == method_cyto(2)/method_cyto(3),...
        'method_cyto(3) must evenly divide method_cyto(2)' );

% Theta-parameter of the theta-method
theta=0.5; % []

% Relaxation parameter in the fixed-point iteration: it must be 0<alpha<=1   (alpha=1 means no relaxation)
alpha=1; % []

% Relative error on the fixed-point iteration of the homogenized problem 
tol_fix=1e-5; % []

% Type of norm used to estimate convergence: (L^\infty (if true), or L^1 (if false))
norma_inf=false;

% Duration of simulation
t_fin = .5; % s

% Number of integration steps
n_step_t = 3*450;

% Downsampling: save the solution at one time step each downsample time steps
%  Useful to reduce the amount of memory required to store sol
%  Reconstruct by linear interpolation.  
%  downsample must divide n_step_t
downsample=1;
n_step_t=downsample*ceil(n_step_t/downsample);

% Record datac
datamat = [datamat;...
         method_cyto(:);...
         theta;...
         alpha;...
         tol_fix;...
         norma_inf;...
         t_fin;...
         n_step_t;...
         downsample];

% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Steady state solver
% maximum unbalancing of fluxes allowed in the steady-state 
tol_stat=1e-8; % [(\mu m) (\mu M)/s = 10^(-9) mole/(m^2 s)]

% Record datac
datamat = [datamat;...
         tol_stat];

% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Physical constants
% Avogadro's number
N_Av = 6.022e23; % #/mol

% Faraday constant
F=96500/1e21; % C/(uM um^3)

% Record datac
datamat = [datamat;...
         N_Av;...
         F];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Cytoplasmic buffers
% cG buffer
B_cG = 1.8036; % Taken from Shen

% Ca buffer
B_Ca = 25.276; % Taken from Shen

% Record datac
datamat = [datamat;...
         B_cG;...
         B_Ca];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% R_st biochemistry
% Activation rate of G_st by R_st
nu_RG   = 211.77; % um^2 /s
               % Note: this number times R_st represents percentage of
               % available unactivated G [#/um^2] that is being turned-over
               % to G* per second

% Diffusion coefficient of R_st
D_R_st = 1.7053; % um^2/s

% Decay rate of R_st
k_R = 87.125; % Hz

% Total surface density of activated and unactivated R
R_sigma = 10e3; % #/um^2

% Record datac
datamat = [datamat;...
         nu_RG;...
         D_R_st;...
         k_R;...
         R_sigma];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% G_st Biochemistry
% Production rate of E_st by G_st
k_GE = .334; % #^(-1)s^(-1) (=s^-1)
          % Note: coefficient is 1 because #G:#PDE is 2:1 meaning #G:#E is
          % 1:1. And so k_GE*G_st = |G_st| *Hz.  It follows the production
          % term for E/1G_st is (E_sigma - E_st)*|G_st|*Hz/|G_st| 
          % This way G_st turns over precisely the E which is available,
          % no more or less.

% Diffusion coefficient of G_st
D_G_st = 1.4172; % um^2/s

% Total surface density of activated and unactivated G_st
G_sigma = 253.82; % #/um^2

% Record datac
datamat = [datamat;...
         k_GE;...
         D_G_st;...
         G_sigma];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% PDE_st biochemistry (free parameters)
% Diffusion coefficient for E_st
D_E_st = 1.47; % um^2/s

% Decay constant of E_st
k_E = 82.331; % Hz

% Surface density of PDE
PDE_sigma = 114.68; % #/um^2

% Dark hydrolysis of cG by E
Beta_dark = 8.5725; % #/s
               % Not used in bright light sims except to compute
               % k_hyd below
               
% Michaelis-Menten constant (for Td:PDE?)
K_m = 10e-6; % M
             % Value from Shen
             
% Turnover rate of single-activated PDE holomer with 1 transducin (=half-value of paper-table's)
k_cat = .33522*1.3351e-5; % (1e-18)mol/s
              % Value so that k_st = 0.9 um^3/s given K_m
              %[?mol]/s/[?] = um^3/s => [?] = [?mol]/s / (um^3/s)
              %[?] = [?mol]/(um^3)
              %If [?] = M => [?] = mol/m^3 = mol/um^3*1e-18
              
% Record datac
datamat = [datamat;...
         D_E_st;...
         k_E;...
         PDE_sigma;...
         Beta_dark;...
         K_m;...
         k_cat];

% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% PDE_st biochemistry (derived parameters)
% Surface density of E
E_sigma = 2*PDE_sigma;

% Volumic density of E
E_vol = E_sigma/(.5*nu*epsilon_0);

% Dark hydrolysis of cG by PDE
k_hyd=0.5*nu*epsilon_0*Beta_dark/PDE_sigma;  % um^3/s
                                             % Actual param in weak
                                             % formulation for bright light

% Hydrolytic efficacy of 1/2-activated PDE dimer (ie E*)
kcat_DIV_Km = k_cat/K_m; % um^3/s;
                         % to convert to uM^-1*s^-1 multiply by 6.022e2
                         

% Rate of hydrolysis of cGMP by light-activated E*
k_st=kcat_DIV_Km/B_cG; % um^3/(s)
                       % Shen has k_st = 0.9 
                       % Formula is (k_cat/K_m)/(B_cG)
                       % Dropped the N_Av shown in Shen bc kcat/K_m was given in 
                       % mol/s not #/s, the latter which N_Av converts
                       
% Record datac
datamat = [datamat;...
         E_sigma;...
         E_vol;...
         k_hyd;...
         kcat_DIV_Km;...
         k_st];

% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% cG biochemistry
% Trial cG dark value
u_tent = 2; % uM

% Diffusion coefficient of cG
kk_u = 89.216; % um^2/s

% Record datac
datamat = [datamat;...
         u_tent;...
         kk_u];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Ca biochemistry
% Trial Ca dark value
v_tent = .5; % uM

% Diffusion coefficient of Ca
kk_v = 15; % um^2/s

% Record datac
datamat = [datamat;...
         v_tent;...
         kk_v];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% cG-gated current biochemistry
% Saturating cG-gated current 
%  Prenormalize by 1/(B_Ca*F) due to roundoff concerns
j_cG_max=3138.2e-12; % [A]
j_cG_max=j_cG_max/(B_Ca*F);


% Hill coefficient
m_cG=2.9228;

% Half-maximum activation concentration, assumed constant
K_cG=14.177; % uM

% fraction of cGMP-activated current given by Ca2+
f_Ca = 0.25577;

% Record datac
datamat = [datamat;...
         j_cG_max;...
         m_cG;...
         K_cG;...
         f_Ca];
    
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Exchanger current biochemistry
% Saturating exchange current
j_ex_sat=8.7861e-12; % [A]
j_ex_sat=j_ex_sat/(B_Ca*F);

% Half-maximum activation concentration
K_ex = 0.40262; % uM

% Record datac
datamat = [datamat;...
         j_ex_sat;...
         K_ex];
    
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% Cyclase biochemistry
% Maximum rate of production of cGMP
alpha_max = 55.802; % uM/s

% Minimum rate of cGMP production
alpha_min = alpha_max/5.5076; % uM/s

% Hill coefficient
m_cyc = 2.3027;

% Half-maximum activation concentration
K_cyc = .13424; % uM

% Record datac
datamat = [datamat;...
         alpha_max;...
         alpha_min;...
         m_cyc;...
         K_cyc];
     
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];
       
%% Illumination Parameters
% Number of uniformly distributed R* molecules at time t = 0
%  Set value to NaN if a space varying concentration is intended and define
%  this function in cyto/cascade/Rst_0.m
n_Rst0 = 940; % Use bleach = eta0*(.7*(p(3) >= 6).*(p(3)<=7)); to get Shen SPR

% Continuous flash intensity rate for steady illumination
%  Set value to NaN if a space or time varying concentration is intended
%  and define this function in cyto/cascade/flash.m
rate = 0; % #VP activated/s if at dark adapted state

% Record into datac
datamat = [datamat;...
           n_Rst0;...
           rate];
       
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

%% ODE Restart Parameters
% Data needed to restart a simulation for a partially simulated system
flag_restart = false;

% The last completed,simulated column index of taxis
rst_index = 0;

% Record into datac
datamat = [datamat;...
           flag_restart;...
           rst_index];
       
% Add into pointer
pointer = [pointer;...
           size(datamat,1)];

% Make character name for restart files
act_fname = 'activation.mat';
Evol_fname= 'Evol.mat';
msg_fname = 'messenger.mat';

fnames = {act_fname,Evol_fname,msg_fname};
end

