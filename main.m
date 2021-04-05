% Batch script for simulating cone bright light response. 
%% Indicate the type of job
% flag_job may be 'Act','Evol','Msg','SA_pd','SA_Sbl','Fit','RanW'
assert(exist('flag_job') == 1,...
       'flag_job must be defined before running main.m');

% Check if M_gl, K_gl have been preassembled
if exist('flag_preassembled') ~= 1
    flag_preassembled = 0;
end

%% Run job in a batch folder
% Make the directory
mkdir batch

% cd into batch
cd batch

%% Copy files into batch and run appropriate job

% Copy needed files into batch
copyfile ..\common\* .\
copyfile ..\cyto\*.m .\
copyfile ..\elements\*.m .\
copyfile '..\ode solver'\* .\


switch flag_job
    %% Activation
    case 'Act'
        copyfile ..\cyto\Cascade\* .\
        copyfile ..\job_act.m .\
        
        % Grab the data
        [datamat,pointer] = data_set;
        
        % Run the job
        job_act
        
        % Copy data back into home directory
        movefile activation.mat ..\
        
        % Delete the batch folder
        cd ..\
        delete batch\*
        rmdir batch\
        
        % clear workspace
        clear all
        
    
    %% Evolume
    case 'Evol'
        copyfile '..\cyto\2nd Messenger\*' .\
        copyfile ..\job_Evol.m .\
        
        % Copy in needed simulation data 
        copyfile '..\activation.mat' .\
        
        % Grab the data 
        [datamat,pointer] = data_set;
        
        % Run the job
        job_Evol
        
        % Copy data back into home directory
        movefile Evol.mat ..\
        
        % Delete the batch folder
        cd ..\
        delete batch\*
        rmdir batch\
        
        % clear workspace
        clear all
        
    %% Messenger
    case 'Msg'
        copyfile '..\cyto\2nd Messenger\*' .\
        copyfile ..\job_msg.m .\
        
        % Copy in needed simulation data
        copyfile ..\Evol.mat .\
        
        % Grab the data
        [datamat,pointer] = data_set;
        
        % Run the job
        job_msg
        
        % Copy data back into home directory
        movefile messenger.mat ..\
        
        % Delete the batch folder
        cd ..\
        delete batch\*
        rmdir batch\
        
        % clear workspace
        clear all
        
    %% Local Sensitivity Analysis
    case 'SA_pd'
        disp(['If restarting a parallel run, the progress.mat file ' newline ...
              'should be in parent directory with main.m']);
        disp(['If job did not complete, you should manually transfer the ' newline ...
              'trials and progress.mat file back to the parent directory']);
          
        % Transfer remaining files into batch 
        copyfile ..\cyto\Cascade\*.m .\
        copyfile ..\cyto\'2nd Messenger'\*.m .\
        copyfile ..\'sensitivity analysis'\*.m .\
        copyfile ..\job_SApd.m .\
        
        % Setup progress.mat file if not already there
        try
           movefile ..\progress.mat .\
        catch errorz
            flag_base = true;
            flag_restart = false;
            save('progress.mat');
        end
        clear errorz
        
        % Run the job
        job_SApd
        
        % Transfer files back to main directory
        movefile .\*.mat ..\
        
        % Delete the batch folder
        cd ..\
        delete batch\*
        rmdir batch\
        
        % Clear the workspace
        clear all
    
    %% Parameter fitting by random sampling
    case 'Fit'
        disp(['Variable n_trials should have been set by user already' ...
              newline ...
              'Variable flag_geomvary should be set by user already to ' ...
              'true/false to say if Rb,Rt,H participate in the fitting']);
        disp(['Note: data_set.m is only used for params that are not ' ...
              'randomly sampled. '...
              newline ...
              'The distributions for others are defined solely in mc_sample.m'...
              newline ...
              'To flag subset of params to vary use 2nd section of job_fit']);
        
        % Transfer remaining files into batch 
        copyfile ..\cyto\Cascade\*.m .\
        copyfile ..\cyto\'2nd Messenger'\*.m .\
        copyfile ..\'sensitivity analysis'\*.m .\
        copyfile ..\'sensitivity analysis'\'param fit'\*.m .\
        copyfile ..\job_fit.m .\
        try 
            copyfile ..\RandNum.mat .\
            flag_rngrst = true;
        catch
            flag_rngrst = false;
        end
        
        % Run the job and set n_trials
        job_fit
        
        % Transfer files back to main directory
        movefile .\*.mat ..\
        
        % Delete the batch folder
        cd ..\
        delete batch\*
        rmdir batch\
        
        % Clear the workspace
        clear all
    
    %% Global sensitivity analysis by Sobol indices
    case 'SA_Sbl'
        disp('prepr_Sbl should have been run already with trial folder here');
        
        % Transfer remaining files into batch 
        copyfile ..\cyto\Cascade\*.m .\
        copyfile ..\cyto\'2nd Messenger'\*.m .\
        copyfile ..\'sensitivity analysis'\*.m .\
        copyfile ..\'sensitivity analysis'\'Sobol'\*.m .\
        copyfile ..\job_SASbl.m .\
        copyfile ..\trials\* .\
        
        % Run the job and set n_trials
        job_SASbl
        
        % Transfer files back to main directory
        movefile .\*.mat ..\
        movefile ..\trial*mat ..\trials\
        
        
        % Delete the batch folder
        cd ..\
        delete batch\*
        rmdir batch\
        
        % Update OSC job status file
        fileID = fopen('jobstat_osc.txt','w');
        fprintf(fileID,'%d',n_complete);
        fclose(fileID);
        
        % Clear the workspace
        clear all
    
    %% Parameter fitting by random walk
    case 'RanW'
        disp(['Variable n_trials should have been set by user already'...
               newline ...
              'Variable flag_geomvary should be set by user already to ' ...
              'true/false to say if Rb,Rt,H participate in the fitting']);
        disp(['Note: data_set.m is only used for params that are not ' ...
              'randomly sampled. '...
              newline ...
              'The distributions for others are defined solely in mcmc_ranw.m'...
              newline ...
              'To flag subset of params to vary use 2nd section of job_ranw'...
              newline ...
              'To give absolute constraint on ranges use 2nd section of mcmc_alpha'...
              newline ...
              'Note wellposed_ss does enforce some additional constraints']);
        
        % Transfer remaining files into batch 
        copyfile ..\cyto\Cascade\*.m .\
        copyfile ..\cyto\'2nd Messenger'\*.m .\
        copyfile ..\'sensitivity analysis'\*.m .\
        copyfile ..\'sensitivity analysis'\'param fit'\*.m .\
        copyfile ..\job_ranw.m .\
        copyfile ..\cyto\Cascade\* .\
        copyfile ..\job_act.m .\
        copyfile '..\cyto\2nd Messenger\*' .\
        copyfile ..\job_Evol.m .\
        copyfile '..\cyto\2nd Messenger\*' .\
        copyfile ..\job_msg.m .\
        try 
            copyfile ..\RandNum.mat .\
            flag_rngrst = true;
        catch
            flag_rngrst = false;
        end
        
        % Run the job and set n_trials
        job_ranw
        
        % Transfer files back to main directory
        movefile .\*.mat ..\
        
        % Delete the batch folder
        cd ..\
        delete batch\*
        rmdir batch\
        
        % Clear the workspace
        clear all
        
    otherwise
        error('flag_job does not match a valid job type')
end



