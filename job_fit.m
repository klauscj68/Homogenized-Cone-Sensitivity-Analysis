% Run a randomized-fitness based fitting scheme of parameters 
%% Generate a datamat and pointer
[datamat,pointer,...
          fnames] = data_set();
      
% Check that base values are well-posed in their steady state
[~,flag_solvable] = wellposed_ss(datamat,pointer);
assert(flag_solvable == true,...
       'The base data file did not have well-posed steady state');
   
%% Pregenerate base M matrices for admissible steady state speed
if flag_geomvary == false
    [R_b,R_t,H,cosgamma0,theta_in,theta_fin,epsilon_0,nu,sigma,...
          flag_ch,...
          n_sez,taglia,tol_R,tol_angle,...
          method_cyto,theta,alpha,tol_fix,norma_inf,t_fin,...
          n_step_t,downsample,...
          tol_stat,...
          N_Av,F,...
          B_cG,B_Ca,...
          nu_RG,D_R_st,k_R,R_sigma,...
          k_GE,D_G_st,G_sigma,...
          D_E_st,k_E,PDE_sigma,Beta_dark,K_m,k_cat,...
          E_sigma,E_vol,k_hyd,kcat_DIV_Km,k_st,...
          u_tent,kk_u,...
          v_tent,kk_v,...
          j_cG_max,m_cG,K_cG,f_Ca,...
          j_ex_sat,K_ex,...
          alpha_max,alpha_min,m_cyc,K_cyc,...
          n_Rst0,rate,...
          flag_restart,rst_index] = ...
          data(datamat,pointer);
    
    [pts,prisms,faces_sl,faces_fo,...
         n_pts,n_prism,n_fsl,n_ffo] = ...
      genmesh(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,sigma,...
              n_sez,taglia,tol_R,tol_angle);
              
    
    M = mvol(pts,prisms,...
             R_b,R_t,H);
    Msl = mbd(pts,faces_sl,...
              R_b,R_t,H);
    Mfo = mbd(pts,faces_fo,...
              R_b,R_t,H);
end
      
%% Specify varied params
flag_params = true(max(size(datamat)),1);

%% Shuffle RNG so ind from any earlier runs
rng('shuffle');
RNgen = rng;

% Save the random number state in case of later replication
save('RandNum.mat','RNgen');

%% Submit the jobs in parallel
k = n_trials;
smp_attempts = 0;

WBar = waitbar(0,'Attempting to sample admissible parameters');
tic
while k >= 1
    % Take a random sample of parameter values
    smp_attempts = smp_attempts + 1;
    [smp] = mc_sample(flag_params,datamat,pointer);
    
    % Check if samples have admissible steady state
    if flag_geomvary == true
        flag_proceed = wellposed_ss(smp,pointer);
    else
        flag_proceed = wellposed_ss(smp,pointer,...
                                    M,Msl,Mfo);
    end
    
    if flag_proceed == true
        % Save smp values to its mat file
        save(['trial_' num2str(k) '.mat'],'smp','pointer');
    
        % Submit each job in parallel
        Par_Run(k) = parfeval(@par_trial,3,...
                             smp,pointer,0,0);
                         
        % Update counter
        k = k-1;
        
        % Update waitbar
        waitbar((n_trials - k)/n_trials,WBar);
    end
end
toc
close(WBar);
disp(['Finished submitting well posed parallel jobs'...
      newline ...
      num2str(n_trials) ' attempts of ' num2str(smp_attempts) ...
      ' were admissible']);

%% Collect jobs
% Save the data sets as they finish
for l=1:n_trials
    % Grab a result from par pool when it finishes
    [k,wk_act,wk_Evol,wk_msg] = fetchNext(Par_Run);
    
    % Save the result in a mat file
    %  filename
    fname = ['trial_' num2str(k) '.mat'];
    
    %  Give vars their proper names
    %   Note: job_SApd script has more complete identity list
    taxis   = wk_act{4};
    
    J_tot = wk_msg{3};
    J_drop = wk_msg{4};
    cG0 = wk_msg{5};
    Ca0 = wk_msg{6};
    
    % Save into file
    load(fname,'smp','pointer');
    save(fname,...
        'smp','pointer',...
        'J_tot','J_drop','cG0','Ca0',...
        'taxis',...
        '-v7.3');
   % Display progress
   if l == floor(l/100)*100
      disp([num2str(l) ' many trials completed...']);
   end
end
