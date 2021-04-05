% Generate all trial data sets needed to perform Sobol analysis
%  This step is performed all at once so that we can be embarassingly
%   parallel when we run jobs on the supercomputer. Doing it by a prestep
%   ensures all the samples came from a single cycle of a random number
%   generator and hence are fairly regarded as independent

%% Print to screen
disp(['Note: Random number generator state is set internally based on '...
      'whether there is a RandNum.mat file']);
disp(['Note: Distributions are defined solely in mc_sample_Sbl.m'...
       newline ...
      'To flag for sensitivity use 4th section of prepr_Sbl']);  

%% Initialize
% Range of trials being performed
trials_id = [1 10e3];

%% Setup batch directory
mkdir batch
cd batch

% Copy needed files into batch
copyfile ..\common\* .\
copyfile ..\cyto\*.m .\
copyfile ..\elements\*.m .\
copyfile '..\ode solver'\* .\

copyfile ..\cyto\Cascade\*.m .\
copyfile ..\cyto\'2nd Messenger'\*.m .\
copyfile ..\'sensitivity analysis'\*.m .\
copyfile ..\'sensitivity analysis'\'Sobol'\*.m .\
copyfile ..\job_SASbl.m .\

[datamat,pointer] = data_set();
%% Specify subset of params whose variation is measured
% RECALL: alpha_min is normalized to independent [0,1] according to max
%         value permitted by other model parameters for dark steady state.
%         This ensures a parameter set of truly independent rv's so that
%         any arbitrarily chosen set may be varied.  Implicit in this
%         assumption is that all variables are uniform over their
%         intervals, since then the conditional probability of alpha_min
%         with respect to other vars is the joint density (unif) over
%         marginal (constant in alpha_min), hence is again uniform and so
%         is that induced by Z uniform on [0,1] times the Amax value of
%         mc_sample (which is the upper limit to permissible values of
%         alpha_min for existence of steady state assuming j_ex_sat and
%         f_Ca/2*j_cG_max are separated)

flag_params = false(max(size(datamat)),1);
n_params = length(flag_params);

% Specify the parameters which are varied
% Geometry
flag_params(1:3) = true; % R_b,R_t,H
flag_params(6:9) = true; % theta_fin,epsilon_0,nu,sigma

% Cytoplasmic buffers
pos = pointer(5);
flag_params(pos+1:pos+2) = true; % B_Ca,B_cG

% R_st biochemistry
pos = pointer(6);
flag_params(pos+1:pos+4) = true; % nu_RG,D_R_st,k_R,R_sigma

% G_st biochemistry
pos = pointer(7);
flag_params(pos+1:pos+3) = true; % k_GE,D_G_st,G_sigma

% PDE_st biochemistry
pos = pointer(8);
flag_params(pos+1:pos+4) = true; % D_E_st,k_E,PDE_sigma,Beta_dark

pos = pointer(9);
flag_params(pos+4) = true; % k_cat/K_m

% cG biochemistry
pos = pointer(10);
flag_params(pos+2) = true; % kk_u

% Ca2+ biochemistry
pos = pointer(11);
flag_params(pos+2) = true; % kk_v

% cG-gated current biochemistry
pos = pointer(12);
flag_params(pos+1:pos+4) = true; % j_cG_max,m_cG,K_cG,f_Ca

% Exchanger current biochemistry
pos = pointer(13);
flag_params(pos+1:pos+2) = true; % j_ex_sat,K_ex

% Cyclase biochemistry
pos = pointer(14);
flag_params(pos+1:pos+4) = true; % alpha_max,alpha_min,m_cyc,K_cyc

% Free parameters
id_free = find(flag_params);
n_frparam = nnz(flag_params);

%% Restart RNG so ind from earlier runs
try 
    copyfile ..\RandNum.mat .\
    flag_rngrst = true;
catch
    flag_rngrst = false;
end

if flag_rngrst == true
    load('RandNum.mat','RNgen');
    rng(RNgen);
else
    rng('shuffle');
    RNgen = rng; %#ok<NASGU>
    RNgen0 = rng;
    
    % Save the random number states in case of later replication
    save('RandNum.mat','RNgen','RNgen0');
end

%% Generate M1,M2 sample matrices
n_trials = trials_id(2) - (trials_id(1) - 1);
shift = trials_id(1) - 1;

% Build M1
k = n_trials;
smp_attempts = 0;

WBar = waitbar(0,'Generating M1 Trials');
tic
while k >= 1
    % Make first sample of parameters
    smp_attempts = smp_attempts + 1;
    smp = mc_sample_Sbl(flag_params,datamat,pointer);
    
    % Check if samples have admissible steady state
    [~,flag_solvable] = wellposed_ss(smp,pointer);
    assert(flag_solvable == true,...
           ['The alpha_min reparameterization didnt lead to solvable '...
            'steady state']);
        
    % Write into M1
    M1(:,k) = smp; %#ok<SAGROW>
                     
    % Update counter
    k = k-1;
    
    % Update waitbar
    waitbar((n_trials - k)/n_trials,WBar);
end
toc
close(WBar);

% Build M2
k = n_trials;
smp_attempts = 0;

WBar = waitbar(0,'Generating M2 Trials');
tic
while k >= 1
    % Make first sample of parameters
    smp_attempts = smp_attempts + 1;
    smp = mc_sample_Sbl(flag_params,datamat,pointer);
    
    % Check if samples have admissible steady state
    [~,flag_solvable] = wellposed_ss(smp,pointer);
    assert(flag_solvable == true,...
           ['The alpha_min reparameterization didnt lead to solvable '...
            'steady state']);
        
    % Write into M1
    M2(:,k) = smp; %#ok<SAGROW>
                     
    % Update counter
    k = k-1;
    
    % Update waitbar
    waitbar((n_trials - k)/n_trials,WBar);
end
toc
close(WBar);

% Save the RN generator state for restarting future runs
RNgen = rng;
save('RandNum.mat','RNgen','-append');

%% Generate the trial mat files for sensitivity analysis like Satelli02

% Columns of M1
WBar = waitbar(0,'Columns of M1');
for i=1:n_trials
    waitbar(i/n_trials,WBar);
    datamat = M1(:,i);
    flag_done = false;
    save(['trial_M1_col' num2str(shift+i) '.mat'],...
          'datamat','pointer','flag_done','-v7.3');
end
close(WBar);

% Columns of M2
WBar = waitbar(0,'Columns of M2');
for i=1:n_trials
    waitbar(i/n_trials,WBar);
    datamat = M2(:,i);
    flag_done = false;
    save(['trial_M2_col' num2str(shift+i) '.mat'],...
          'datamat','pointer','flag_done','-v7.3');
end
close(WBar);

% Columns of Ai
WBar = waitbar(0,'Cycling through Ai');
for i=1:n_frparam
    waitbar(i/n_frparam);
    % Build Ai
    ram = M2;
    ram(id_free(i),:) = M1(id_free(i),:);
    
    % Generate columns of Ai
    for j=1:n_trials
        datamat = ram(:,j);
        flag_done = false;
        save(['trial_A' num2str(i) '_col' num2str(shift+j) '.mat'],...
              'datamat','pointer','flag_done','-v7.3');
    end
    
end
close(WBar)

% Columns of Ani
WBar = waitbar(0,'Cycling through Ani');
for i=1:n_frparam
    % Build Ani
    waitbar(i/n_frparam);
    ram = M1;
    ram(id_free(i),:) = M2(id_free(i),:);
    
    % Generate columns of Ani
    for j=1:n_trials
        datamat = ram(:,j);
        flag_done = false;
        save(['trial_An' num2str(i) '_col' num2str(shift+j) '.mat'],...
              'datamat','pointer','flag_done','-v7.3');
    end
    
end
close(WBar)

%% Save master data set
% Save master copy of M1,M2
save('trial_master.mat',...
     'M1','M2','n_trials','trials_id',...
     'flag_params','n_params','id_free','n_frparam',...
     '-v7.3');
 
clear all %#ok<CLALL>

%% Transfer trials out of batch
cd ..\
mkdir trials
movefile batch\trial*mat trials\ 
movefile batch\RandNum.mat .\
delete batch\*
rmdir batch\