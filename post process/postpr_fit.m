%Post-process randomized sims by computing their fitness
%% Define path for auxilliary files
oldpath = addpath('..\sensitivity analysis\param fit\');
addpath('..\Simulation Sets\Param_Fit\');
addpath('..\common\');

%% Initialize variables
n_trials = input('Number of trials: ');
gen_fitness = inf(n_trials,1);

% Store flags to say whether trial looks blatently unconvergent
flag_conv = true(n_trials,1);

%% Compute trial fitness
WBar = waitbar(0,'Computing Fitness of Trials');
for k=1:n_trials
    % job_fit creates all trial mats even when they haven't completed, so
    % must code for the cases that didn't finish
    try
        gen_fitness(k) = fitness(k);
        
        % Do a check to see if trial looks blatently unconvergent
        load(['trial_' num2str(k) '.mat'],'J_drop');
        if ( min(J_drop) <= -75 )||( imag(min(J_drop)) >= 1e-5 )
            flag_conv(k) = false;
        end   
        
    catch ME
        true;
    end
    waitbar(k/n_trials,WBar);
end
close(WBar);

%% Sort trials from most fit to least
[sort_fitness,order_fitness] = sort(gen_fitness);

% Reorder flag_conv to match the sorting
flag_conv = flag_conv(order_fitness);

%% Print and plot summary statistics
% Find the valid simulations
flag_complete = ~isinf(sort_fitness);
flag_real = ( imag(sort_fitness) <= 1e-5 );
flag_valid = logical((flag_complete.*flag_real).*flag_conv);
n_complete = nnz(flag_complete);
n_valid = nnz(flag_valid);

% Redefine sort_fitness and order_fitness to only be valid ones
sort_fitness = sort_fitness(flag_valid);
order_fitness = order_fitness(flag_valid);

% Print and Plot
disp([ num2str(n_complete) ' trials of ' num2str(n_trials) ' completed']);
disp([ num2str(n_valid) ' trials of ' num2str(n_trials) ' were valid']);
disp(['Most successful trial: #' num2str(order_fitness(1))]);

H = histogram(real(sort_fitness),'Normalization','probability');

%% Reset path
path(oldpath);



