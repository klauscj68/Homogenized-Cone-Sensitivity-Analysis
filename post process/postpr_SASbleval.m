% Evaluate all functional over all columns of M1,M2,Ai,Ani
% Path on supercomputer for trials
oldpath = addpath('../trials/');
%% Initialize
% Load trial_master for size and datamat values
load('trial_master.mat');
n_functionals = 16;

% Define time ranges for functionals
tact = [0 .01]; % Remember to visualize in exp flash response, these ranges
trec = [.135 .5];% get a +.005 shift

%% Loop over jobids and eval functionals
% Initialize cell containers
ai = cell(n_frparam,1);
ani = cell(n_frparam,1);

for i=1:nbatch
    fname = ['jobid_batch' num2str(i) '.mat'];
    load(fname);
    
    for j=1:length(jobid)
        id = jobid(j);
        
        % Load datamats
        if id == Inf
            % M1
            DATAMAT = M1;
        elseif id == -Inf
            % M2
            DATAMAT = M2;
        elseif id > 0
            % Ai
            DATAMAT = Ai{id};
        elseif id < 0
            % Ani
            DATAMAT = Ani{-id};
        end
        
        ram = zeros(n_functionals,n_trials);
        parSbl_results = Sbl_results{j};
        parfor k=1:n_trials
            R = parSbl_results{k};
            ram(:,k) = postpr_functional(tact,trec,R.taxis,...
                                           R.J_tot,R.J_drop,...
                                           R.mass_Evol,...
                                           DATAMAT(:,k),pointer,...
                                           true);
        end
        
        % Save the result
        if id == Inf
            % M1
            m1 = ram;
        elseif id == -Inf
            % M2
            m2 = ram;
        elseif id > 0
            % Ai
            ai{id} = ram;
        elseif id < 0
            % Ani
            ani{-id} = ram;
        end
        disp(['Finished task ' num2str(j) '/' num2str(length(jobid)) ...
              ' at worker ' num2str(i)]);
    end
end

%% Save to mat file
save('Sbl_Evals','m1','m2','ai','ani','-v7.3');