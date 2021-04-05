% Run the randomized sample of parameters for Sobol analysis
%% Initialize
assert(isfile('trial_master.mat'),...
       'master.mat file must be in directory');
load('trial_master.mat','n_trials','trials_id',...
     'flag_params','n_params','id_free','n_frparam');


%% Generate a datamat and pointer
[~,pointer,...
          ~] = data_set();

%% Submit the jobs in parallel
% For counting progress rate at start of simulation
last = 0;
n_complete = 0;

clock = 0:10:100;
pos = 1;

% Track the shift for trial id's
shift = trials_id(1)-1;

%% Submit M1
for k=n_trials:-1:1
    fname = ['trial_M1_col' num2str(shift+k) '.mat'];
    
    % If file isn't in directory or if it has already been completed cycle
    %  to next value
    if ~isfile(fname)
        continue
    end
    
    load(fname);
    if flag_done
        n_complete = n_complete+1;
        continue
    end
    
    % Run job
    [wk_act,wk_Evol,wk_msg] = par_trial(datamat,pointer,0,0);     %#ok<UNRCH>
                             
    % Save the result in its mat file
    
    %  Give vars their proper names
    %   Note: job_SApd script has more complete identity list
    taxis   = wk_act{4};
    
    mass_Evol = wk_Evol{3};
    
    J_tot = wk_msg{3};
    J_drop = wk_msg{4};
    cG0 = wk_msg{5};
    Ca0 = wk_msg{6};
    
    % Append the variables to the file
    flag_done = true;
    save(fname,...
         'taxis','mass_Evol','J_tot','J_drop','cG0','Ca0','flag_done',...
         '-append');
    if (k/n_trials*100 >= clock(pos+1))
        pos = pos+1;
        disp(['M1 at ' num2str(clock(pos)) '%'])
    end
     
end

%% Submit M2
pos = 1;
for k=n_trials:-1:1
    fname = ['trial_M2_col' num2str(shift+k) '.mat'];
    
    % If file isn't in directory or if it has already been completed cycle
    %  to next value
    if ~isfile(fname)
        continue
    end
    
    load(fname);
    if flag_done
        n_complete = n_complete+1;
        continue
    end
    
    % Run job
    [wk_act,wk_Evol,wk_msg] = par_trial(datamat,pointer,0,0);    
                             
    % Save the result in its mat file
    
    %  Give vars their proper names
    %   Note: job_SApd script has more complete identity list
    taxis   = wk_act{4};
    
    mass_Evol = wk_Evol{3};
    
    J_tot = wk_msg{3};
    J_drop = wk_msg{4};
    cG0 = wk_msg{5};
    Ca0 = wk_msg{6};
    
    % Append the variables to the file
    flag_done = true;
    save(fname,...
         'taxis','mass_Evol','J_tot','J_drop','cG0','Ca0','flag_done',...
         '-append');
     
     if (k/n_trials*100 >= clock(pos+1))
        pos = pos+1;
        disp(['M2 at ' num2str(clock(pos)) '%'])
    end
                   
end

%% Submit Ai
for i=1:n_frparam
    pos = 1;
    for k = n_trials:-1:1
        fname = ['trial_A' num2str(i) '_col' num2str(shift+k) '.mat'];
    
        % If file isn't in directory or if it has already been completed 
        %  cycle to next value
        if ~isfile(fname)
            continue
        end
    
        load(fname);
        if flag_done
            n_complete = n_complete+1;
            continue
        end
    
       % Run job
        [wk_act,wk_Evol,wk_msg] = par_trial(datamat,pointer,0,0);   
                             
        % Save the result in its mat file
    
        %  Give vars their proper names
        %   Note: job_SApd script has more complete identity list
        taxis   = wk_act{4};
    
        mass_Evol = wk_Evol{3};
    
        J_tot = wk_msg{3};
        J_drop = wk_msg{4};
        cG0 = wk_msg{5};
        Ca0 = wk_msg{6};
    
        % Append the variables to the file
        flag_done = true;
        save(fname,...
             'taxis','mass_Evol','J_tot','J_drop','cG0','Ca0','flag_done',...
             '-append');
   
        if (k/n_trials*100 >= clock(pos+1))
            pos = pos+1;
            disp(['A' num2str(i) ' at ' num2str(clock(pos)) '%'])
        end
     end
end

%% Submit Ani

for i=1:n_frparam
    pos = 1;
    for k = n_trials:-1:1
        fname = ['trial_An' num2str(i) '_col' num2str(shift+k) '.mat'];
    
        % If file isn't in directory or if it has already been completed 
        %  cycle to next value
        if ~isfile(fname)
            continue
        end
    
        load(fname);
        if flag_done
            n_complete = n_complete+1;
            continue
        end
    
        % Run job
        [wk_act,wk_Evol,wk_msg] = par_trial(datamat,pointer,0,0);
                             
        % Save the result in its mat file
    
        %  Give vars their proper names
        %   Note: job_SApd script has more complete identity list
        taxis   = wk_act{4};
    
        mass_Evol = wk_Evol{3};
    
        J_tot = wk_msg{3};
        J_drop = wk_msg{4};
        cG0 = wk_msg{5};
        Ca0 = wk_msg{6};
    
        % Append the variables to the file
        flag_done = true;
        save(fname,...
            'taxis','mass_Evol','J_tot','J_drop','cG0','Ca0','flag_done',...
            '-append');
   
        if (k/n_trials*100 >= clock(pos+1))
            pos = pos+1;
            disp(['An' num2str(i) ' at ' num2str(clock(pos)) '%'])
        end
     end
end