% Evaluate all functional over all columns of M1,M2,Ai,Ani
% Path on supercomputer for trials
oldpath = addpath('../../trials/');
%% Initialize
n_frparam = 33;
n_trials = 50e3;
n_functionals = 16;

% Define time ranges for functionals
tact = [0 .01]; % Remember to visualize in exp flash response, these ranges
trec = [.135 .5];% get a +.005 shift

%% M1 evals
tic
for j=n_trials:-1:1
    load(['trial_M1_col' num2str(j) '.mat']);
    m1(:,j) = postpr_functional(tact,trec,taxis,...
                                           J_tot,J_drop,...
                                           mass_Evol,...
                                           datamat,pointer,...
                                           true);
end
save('Sbl_Evals','m1','-v7.3');
disp('Finished m1');
toc
      
%% M2 evals
tic
for j=n_trials:-1:1
    load(['trial_M2_col' num2str(j) '.mat']);
    m2(:,j) = postpr_functional(tact,trec,taxis,...
                                           J_tot,J_drop,...
                                           mass_Evol,...
                                           datamat,pointer,...
                                           true);
end
save('Sbl_Evals','m2','-append','-v7.3');
disp('Finished m2');
toc
      
%% Ai evals
ai = cell(n_frparam,1);
for i=n_frparam:-1:1
    for j=n_trials:-1:1
        load(['trial_A' num2str(i) '_col' num2str(j) '.mat']);
        ram(:,j) = postpr_functional(tact,trec,taxis,...
                                           J_tot,J_drop,...
                                           mass_Evol,...
                                           datamat,pointer,...
                                           true);
    end
    ai{i,1} = ram;
    save('Sbl_Evals','ai','-append','-v7.3');
    disp(['Finished a' num2str(i)]);
    
end

%% Ani evals
ani = cell(n_frparam,1);
for i=n_frparam:-1:1
    for j=n_trials:-1:1
        load(['trial_An' num2str(i) '_col' num2str(j) '.mat']);
        ram(:,j) = postpr_functional(tact,trec,taxis,...
                                           J_tot,J_drop,...
                                           mass_Evol,...
                                           datamat,pointer,...
                                           true);
    end
    ani{i,1} = ram;
    save('Sbl_Evals','ani','-append','-v7.3');
    disp(['Finished an' num2str(i)]);
    
end