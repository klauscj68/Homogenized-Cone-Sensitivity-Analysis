%Post-process sims by computing changes in all functionals
%  If at some point the number of functionals is changed, you will need to
%  adjust the 'Compute the functionals for base parameter set' and the 
%  'Compute the functional values for all the parameter variations' in
%  addition to changes already made in postpr_functional.
%  'Compute Sensitivity Indices' section has a manual .05 factor that needs
%  to match the job_SApd section Setup parfeval var factor
%% Define the time intervals used for measuring activation and recovery
tact = [0 .01];
trec = [.135 .5];

%% Define path for auxilliary files
oldpath = addpath('..\common\');
addpath('..\');
addpath('..\sensitivity analysis\');
addpath('..\cyto\2nd Messenger\');
addpath('..\elements\');

% Define path with data
addpath('..\');

%% Setup A Cell with the Parameter Names and Initialize Sensitivity Storage
load('progress.mat','flag_params_orig','names');

% Functionals are Iact, Ipeak, Iwidth, Estact, Estpeak, Estrec
n_data = size(flag_params_orig,2);
index_params = find(flag_params_orig);
n_params = size(index_params,1)/2;

%% Compute the functionals for base parameter set
% Load data
load('cyto_0.mat')

% Compute all functionals
[datamat,pointer,...
                fnames] = data_set();
fct0 = postpr_functional(tact,trec,taxis,...
                         J_tot,J_drop,...
                         mass_Evol,...
                         datamat,pointer,true);

%% Compute the functional values for all the parameter variations
% Functional value array
%  Each row is a functional of postpr_functional
%  In a single row, the 1:n_data block is the first variation of job_SApd
%  of the Setup parfeval section and n_data + 1:2*n_data is the second
%  variation. That is each column corresponds to a +- var in the matching
%  parameter value
FCT = zeros(16,2*n_data);
WBar = waitbar(0,'Running Parameter Sets');
for k=1:2*n_params
        [j,i] = ind2sub(size(flag_params_orig),index_params(k));
        load(['cyto_' names{i} '_var' num2str(j) '.mat'])
        
        % Compute the new functional values for param i var j
        fct = postpr_functional(tact,trec,taxis,...
                         J_tot,J_drop,...
                         mass_Evol,...
                         smp,pointer,true);
                     
        % Store it in the array
        FCT(:,i+(j-1)*n_data) = fct;
        
    waitbar(k/(2*n_params),WBar);
end
close(WBar)

%% Trim Functionals down to relevant parameters
flag_params = logical(flag_params_orig(1,:));

% Perform trimming
%  Use two copies because FCT has a +- variation in each row
FCT = FCT(:,repmat(flag_params,1,2));

%% Compute the Sensitivity Indices for the Functionals
% Grab stepsizes used in parameter increment from base data sheet
[datamat,pointer,...
          fnames] = data_set();

prmtbase = datamat(flag_params);
prmtbase = prmtbase'; % Switch to row

prmtstep = 1./(.05*prmtbase);

% Compute the partials
%  The 1:n_params block of columns is bdiff
%  The n_params+1:2*n_params block of columns is fdiff
PD = ([-prmtstep prmtstep].*(FCT - fct0));

%  Average forward and backwards difference quotients
%   Avg encodes the weights used on resp bdiff and fdiff
Avg = [0 1];
PD = Avg(1)*PD(:,1:n_params) + Avg(2)*PD(:,n_params + (1:n_params));

% Normalize for relative sensitivity: (dy/dx)/(y/x)
%  Careful, bc MATLAB implicity expands the singleton dim by repmat ./* is
%   not associative here as this repmat happens in different dimensions for
%   fct0 (col vector) and prtmbase (row vector)
%  Like PD, this resembles a Jacobian matrix
SI = (PD./fct0).*prmtbase;

%% Declare functional names for table
% Transpose is for putting in table format
SIact   =  SI(1,:)';
SIpeak  =  SI(2,:)';
SIrec   =  SI(3,:)';
SEstact =  SI(4,:)';
SEstpeak =  SI(5,:)';
SEstrec =  SI(6,:)';
SL2 =  SI(7,:)';
STpeak = SI(8,:)';
SIwidth = SI(9,:)';
SJdark = SI(10,:)';
SJovsht = SI(11,:)';

%% Construct a table with the SI indices
SI_Table = table(SIact,SIpeak,SIrec,SEstact,SEstpeak,SEstrec,SJdark,...
                 SJovsht,SL2,STpeak,SIwidth);
ROWNAMES = names(flag_params);
ROWNAMES = reshape(ROWNAMES,n_params,1);
SI_Table.Rows = ROWNAMES;
SI_Table = [SI_Table(:,12) SI_Table(:,1:11)];


%% Save the table
writetable(SI_Table,'SI_pd.csv')

%% Reset path
path(oldpath);