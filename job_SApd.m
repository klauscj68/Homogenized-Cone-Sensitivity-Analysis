% Run parameter sensitivity analysis on bright light cone model
disp(['Before beginning you should comment out ode_integrate.m and steady_state.m' newline ...
      ' printing time integration progress and steady states to std output.' newline ...
      'Also comment out ode_integrate.m saving the global matrices']);
disp(['There should be a progress.mat file with a flag_restart set to ' newline ...
      'true or false and if true a flag_params value there to overwrite below']);
disp(['progress.mat should also contain a flag_base set to true or false ' newline ...
      'depending on whether the base values of parameters should be run']);

%% Define path for auxilliary files
% No longer using this section because assuming user will run trials using
% main.m which handles files in different way

%oldpath = addpath('common\');
%addpath('cyto\');
%addpath('cyto\Cascade\');
%addpath('cyto\2nd Messenger\');
%addpath('elements\');
%addpath('ode solver\');
%addpath('sensitivity analysis\');

%% Grab the base data
[datamat,pointer] = data_set;
  
%% Set logical flag for params wish to measure sensitivity
% Define character cell with the parameter names for saving data
n_flag = size(pointer,1);
names = cell(1,n_flag);
flag_params = zeros(1,n_flag);

names{1} = 'R_b'; flag_params(1) = 1;
names{2} = 'R_t'; flag_params(2) = 1;
names{3} = 'H'; flag_params(3) = 1;
names{4} = 'cosgamma0'; flag_params(4) = 0;
names{5} = 'theta_in';  flag_params(5) = 0;
names{6} = 'theta_fin'; flag_params(6) = 1;
names{7} = 'epsilon_0'; flag_params(7) = 1;
names{8} = 'nu'; flag_params(8) = 1;
names{9} = 'sigma'; flag_params(9) = 1;
names{10} = 'flag_ch'; flag_params(10) = 0;
names{11} = 'flag_ch'; flag_params(11) = 0;
names{12} = 'flag_ch'; flag_params(12) = 0;

B = pointer(1); % Last block completed
names{B + 1} = 'n_sez'; flag_params(B + 1) = 0;
names{B + 2} = 'taglia'; flag_params(B + 2) = 1;
names{B + 3} = 'tol_R'; flag_params(B + 3) = 1;
names{B + 4} = 'tol_angle'; flag_params(B + 4) = 0;

B = pointer(2);
names{B + 1} = 'method_cyto'; flag_params(B + 1) = 0;
names{B + 2} = 'method_cyto'; flag_params(B + 2) = 0;
names{B + 3} = 'method_cyto'; flag_params(B + 3) = 0;
names{B + 4} = 'theta'; flag_params(B + 4) = 0;
names{B + 5} = 'alpha'; flag_params(B + 5) = 0;
names{B + 6} = 'tol_fix'; flag_params(B + 6) = 0;
names{B + 7} = 'norma_inf'; flag_params(B + 7) = 0;
names{B + 8} = 't_fin'; flag_params(B + 8) = 1;
names{B + 9} = 'n_step_t'; flag_params(B + 9) = 0;
names{B + 10} = 'downsample'; flag_params(B + 10) = 0;

B = pointer(3);
names{B + 1} = 'tol_stat'; flag_params(B + 1) = 1;

B = pointer(4);
names{B + 1} = 'N_Av'; flag_params(B + 1) = 0;
names{B + 2} = 'F'; flag_params(B + 2) = 0;

B = pointer(5);
names{B + 1} = 'B_cG'; flag_params(B + 1) = 1;
names{B + 2} = 'B_Ca'; flag_params(B + 2) = 1;

B = pointer(6);
names{B + 1} = 'nu_RG'; flag_params(B + 1) = 1;
names{B + 2} = 'D_R_st'; flag_params(B + 2) = 1;
names{B + 3} = 'k_R'; flag_params(B + 3) = 1;
names{B + 4} = 'R_sigma'; flag_params(B + 4) = 1;

B = pointer(7);
names{B + 1} = 'k_GE'; flag_params(B + 1) = 1;
names{B + 2} = 'D_G_st'; flag_params(B + 2) = 1;
names{B + 3} = 'G_sigma'; flag_params(B + 3) = 1;

B = pointer(8);
names{B + 1} = 'D_E_st'; flag_params(B + 1) = 1;
names{B + 2} = 'k_E'; flag_params(B + 2) = 1;
names{B + 3} = 'PDE_sigma'; flag_params(B + 3) = 1;
names{B + 4} = 'Beta_dark'; flag_params(B + 4) = 1;
names{B + 5} = 'K_m'; flag_params(B + 5) = 1;
names{B + 6} = 'k_cat'; flag_params(B + 6) = 1;

B = pointer(9);
names{B + 1} = 'E_sigma'; flag_params(B + 1) = 1;
names{B + 2} = 'E_vol'; flag_params(B + 2) = 1;
names{B + 3} = 'k_hyd'; flag_params(B + 3) = 1;
names{B + 4} = 'kcat_DIV_Km'; flag_params(B + 4) = 1;
names{B + 5} = 'k_st'; flag_params(B + 5) = 1;

B = pointer(10);
names{B + 1} = 'u_tent'; flag_params(B + 1) = 1;
names{B + 2} = 'kk_u'; flag_params(B + 2) = 1;

B = pointer(11);
names{B + 1} = 'v_tent'; flag_params(B + 1) = 1;
names{B + 2} = 'kk_v'; flag_params(B + 2) = 1;

B = pointer(12);
names{B + 1} = 'j_cG_max'; flag_params(B + 1) = 1;
names{B + 2} = 'm_cG'; flag_params(B + 2) = 1;
names{B + 3} = 'K_cG'; flag_params(B + 3) = 1;
names{B + 4} = 'f_Ca'; flag_params(B + 4) = 1;

B = pointer(13);
names{B + 1} = 'j_ex_sat'; flag_params(B + 1) = 1;
names{B + 2} = 'K_ex'; flag_params(B + 2) = 1;

B = pointer(14);
names{B + 1} = 'alpha_max'; flag_params(B + 1) = 1;
names{B + 2} = 'alpha_min'; flag_params(B + 2) = 1;
names{B + 3} = 'm_cyc'; flag_params(B + 3) = 1;
names{B + 4} = 'K_cyc'; flag_params(B + 4) = 1;

B = pointer(15);
names{B + 1} = 'n_Rst0'; flag_params(B + 1) = 0;
names{B + 2} = 'rate'; flag_params(B + 2) = 0;

B = pointer(16);
names{B + 1} = 'flag_restart'; flag_params(B + 1) = 0;
names{B + 2} = 'rst_index'; flag_params(B + 2) = 0;

assert(size(datamat,1) == B + 2,'Some parameters in datasheet were not named');
               
% Duplicate flags for two variations we want
flag_params = repmat(flag_params,2,1);
               
% Make an original version of flag_params to be called in post processing
flag_params_orig = flag_params;

% Do check to make sure it is formatted like postpr_sims expects
check = (flag_params_orig(1,:) == flag_params_orig(2,:));
check = prod(check);
assert(logical(check),...
       'postpr_sims.m expects forwards and backwards quotients');

%% Specify any already completed parameters
load('progress.mat','flag_restart','flag_base');
if flag_restart == true
    load('progress.mat','flag_params');
end

% Check if want to run the non-varied parameter set
if flag_base == true
    smp = datamat;
   [wk_act,wk_Evol,wk_msg] = par_trial(smp,pointer,...
                                       4,0);
   
   % Save the result in a mat file
    fname = 'cyto_0.mat';
    
    %  Give vars their proper names
    Rst_sig = wk_act{1};
    Gst_sig = wk_act{2};
    Est_sig = wk_act{3};
    taxis   = wk_act{4};
    
    Evol = wk_Evol{1};
    mass_Evol = wk_Evol{3};
    
    cG = wk_msg{1};
    Ca = wk_msg{2};
    J_tot = wk_msg{3};
    J_drop = wk_msg{4};
    cG0 = wk_msg{5};
    Ca0 = wk_msg{6};
    
    % Save into file
    save(fname,...
        'Rst_sig','Gst_sig','Est_sig','Evol','mass_Evol',...
        'cG','Ca','J_tot','J_drop','cG0','Ca0',...
        'taxis',...
        'smp','pointer',...
        '-v7.3');
    
    % Update the progress.mat to reflect it is not needed
    flag_base = false;
    save('progress.mat','flag_base','flag_restart');
                          
end
    
%% Setup parfeval
% Index the parameters to be compatible with parfeval
index_params = find(flag_params);
n_params = size(index_params,1);
     
% Submit parallel calls for parameter sensitivity
var = .05*[-1 1];
for k=1:n_params
    [i,j] = ind2sub(size(flag_params),index_params(k));
    delta_r = var(i);
    
    % Submit each job in parallel
    Par_Run(k) = parfeval(@par_trial,3,...
                          datamat,pointer,j,delta_r);
end

% Save the data sets as they finish
WB = waitbar(0,'Collecting Parallel Runs');
for l=1:n_params
    % Grab a result from par pool when it finishes
    [k,wk_act,wk_Evol,wk_msg] = fetchNext(Par_Run);
    
    % Identify which result this is
    [i,j] = ind2sub(size(flag_params),index_params(k));
    
    % Save the result in a mat file
    %  Names say the parameter being changed, and vari says var(i) was the
    %  component being taken at that moment
    fname = ['cyto_' names{j} '_var' num2str(i)  '.mat'];
    
    %  Save corresponding data set
    smp = datamat;
    delta_r = var(i);
    
    smp(j) = (1+delta_r)*smp(j);
    
    %  Give vars their proper names
    Rst_sig = wk_act{1};
    Gst_sig = wk_act{2};
    Est_sig = wk_act{3};
    taxis   = wk_act{4};
    
    Evol = wk_Evol{1};
    mass_Evol = wk_Evol{3};
    
    cG = wk_msg{1};
    Ca = wk_msg{2};
    J_tot = wk_msg{3};
    J_drop = wk_msg{4};
    cG0 = wk_msg{5};
    Ca0 = wk_msg{6};
    
    % Save into file
    save(fname,...
        'Rst_sig','Gst_sig','Est_sig','Evol','mass_Evol',...
        'cG','Ca','J_tot','J_drop','cG0','Ca0',...
        'taxis',...
        'smp','pointer',...
        '-v7.3');
         
    % Update the flag_params to reflect this new trial has been completed
    flag_params(i,j) = 0;
    flag_restart = true;
    save('progress.mat','flag_params','flag_restart','flag_base',...
         'flag_params_orig','names');
     
    
    % Update the waitbar
    waitbar(l/n_params,WB);
    
end

% Close waitbar
close(WB);

%% Reset path
% No longer needed because running script through main.m
%path(oldpath);

