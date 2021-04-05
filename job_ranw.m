% Run a random walk with fitness scoring scheme parameter fitting
%% Generate a datamat and pointer
[datamat,pointer,...
          fnames] = data_set();
      
%% Specify varied params
flag_params = true(max(size(datamat)),1);

%{ 
% Parameters are varying
B = pointer(6);
flag_params(B+3) = true; % k_R

B = pointer(8);
flag_params(B+2) = true; % k_E
flag_params(B+4) = true; % Beta_dark

B = pointer(9);
flag_params(B+4) = true; % kcat/Km

B = pointer(12);
flag_params(B+1) = true; % j_cG_max
flag_params(B+2) = true; % m_cG
flag_params(B+3) = true; % K_cG

B = pointer(13);
flag_params(B+1) = true; % j_ex_sat
flag_params(B+2) = true; % K_ex

B = pointer(14);
flag_params(B+1) = true; % amax
flag_params(B+2) = true; % amax/amin
flag_params(B+3) = true; % m_cyc
flag_params(B+4) = true; % K_cyc
%}

% Free parameters not varying
flag_params(1:3) = false; % R_b,R_t,H
flag_params(6:9) = false; % theta_fin,epsilon_0,nu,sigma

B = pointer(1);
flag_params(B+1:B+2) = false; % n_sez,taglia

B = pointer(2);
flag_params(B+8:B+9) = false; % t_fin,n_step_t


B = pointer(11);
flag_params(B+2) = false; % kk_v

% Don't vary determined parameters
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

%% Restore RNG so ind from any earlier runs
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

%% Initialize and run base set
n_smp = max(size(datamat));
% To store number of times chain stayed fixed at value thru rejection
nstay = ones(n_trials+1,1); 
% To store parameter values at this position
wsmp = zeros(n_smp,n_trials+1);
wsmp(:,1) = datamat;

% To record quality of sample 
fitscore = NaN(1,n_trials);

% Simulate activation
job_act;
    
% Simulate Evol
job_Evol;
    
% Simulate Msg
job_msg;
    
% Compute fitness 
movefile('messenger.mat',['trial_' num2str(1) '.mat']);
smp = datamat;
save(['trial_' num2str(1) '.mat'],'cG0','Ca0','smp','pointer','-append');

fitscore(1) = fitness(1);

%% Run random walk
tic
%WBar = waitbar(0,'Running random walk');
for k=2:n_trials+1
    flag_accept = false;
    while ~flag_accept
        % Propose the new sample
        flag_propose = false;
        while ~flag_propose
            [propose,prange,threshold,p,pr_gamma] =...
                     mcmc_ranw(flag_params,wsmp(:,k-1),pointer,fitscore(k-1));
        
            % Check if admissible
            if flag_geomvary == false
                flag_propose = wellposed_ss(propose,pointer,...
                                            M,Msl,Mfo);
            else
                flag_propose = wellposed_ss(propose,pointer);
            end
            if ~flag_propose
                % If rejected record so may reconstruct ergodic ensemble
                nstay(k-1) = nstay(k-1)+1;
            end
        end
        datamat = propose;
    
        % Simulate activation
        job_act;
    
        % Simulate Evol
        job_Evol;
    
        % Simulate Msg
        job_msg;
    
        % Compute fitness
        movefile('messenger.mat',['trial_' num2str(k) '.mat']);
        smp = datamat;
        save(['trial_' num2str(k) '.mat'],'cG0','Ca0','smp','pointer','-append');
        newscore = fitness(k);
        
        %% Check if accept
        [flag_accept,flag_nstay] = mcmc_alpha(fitscore(k-1),newscore,...
                                    wsmp(:,k-1),propose,pointer,...
                                    prange,flag_params,...
                                    threshold,p,pr_gamma);
        if flag_accept == true
            % If accepted record the new value
            wsmp(:,k) = propose;
            fitscore(k) = newscore;
        else
            % If rejected record so may reconstruct ergodic ensemble
            nstay(k-1) = nstay(k-1)+1;
        end
    end
    if mod(k,25) == 0
        save('MCMC_Walk.mat','wsmp','fitscore','k','nstay');
        save('RandNum.mat','RNgen','-append');
    end
%    waitbar((k-1)/n_trials,WBar);
end
%close(WBar);
% Save the RN generator state for restarting future runs
RNgen = rng;
save('RandNum.mat','RNgen','-append');
%% Save completed walk and print stats
[best,idx] = min(fitscore);
save('MCMC_Walk.mat','wsmp','fitscore','best','idx','nstay');

%histogram(real(fitscore));
disp(['Best fitness seen is ' num2str(best) ' at trial '...
      num2str(idx)]);
