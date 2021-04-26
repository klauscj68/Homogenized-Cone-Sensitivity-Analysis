% Run the randomized sample of parameters for Sobol analysis
%% Initialize
assert(isfile('trial_master.mat'),...
       'master.mat file must be in directory');
load('trial_master.mat');


%% Generate a datamat and pointer
[~,pointer,...
          ~] = data_set();

%% Submit the jobs in parallel
% For counting progress rate at start of simulation
clock = 0:10:100;
pos = 1;

%% Find worker's load
for i=1:nbatch
    fname = ['jobid_batch' num2str(i) '.mat'];
    if ~isfile(fname)
        continue
    else
        load(fname);
        
        % Loop over mats this worker responsible for
        for k=1:length(jobid)
            j = jobid(k);
            % jobids meaning
            %  +1-n_frparam means Ai{i}
            %  -1-n_frparam means Ani{i}
            %  +Inf means M1
            %  -Inf means M2
            %  NaN means a placeholder that will be skipped
            if isnan(j)||flag_done(k)
                continue
            elseif j == Inf
                ram = M1;
            elseif j == -Inf
                ram = M2;
            elseif j > 0
                ram = Ai{j};
            elseif j < 0
                ram = Ani{-j};
            else
                error('jobids index was not correct');
            end
            
            % Loop over columns of this mat
            if isempty(flag_colpos{k})
                wkcols = 1:n_trials;
                flag_colpos{k} = false(n_trials,1);
            else
                wkcols = 1:n_trials;
                wkcols(flag_colpos{k}) = [];
            end
            
            ncols = length(wkcols);
            % Organize wkcols into small parallel batch sizes for memory
            nloads = 300;
            nrdcols = ceil(ncols/nloads)*nloads;
            wkbatch = [wkcols NaN(1,nrdcols - ncols)];
            wkbatch = reshape(wkbatch,nloads,nrdcols/nloads);
            for l=wkbatch
                for m=nnz(~isnan(l)):-1:1
                % Submit each job in parallel
                Par_Run(m) = parfeval(@par_trial,3,...
                                      ram(:,l(m)),pointer,NaN,0);
                end
                
                for m=1:nnz(~isnan(l))
                    % Fetch jobs and save partial progress
                    [idm,wk_act,wk_Evol,wk_msg] = fetchNext(Par_Run);
                
                    %  Note: job_SApd script has more complete identity list
                    nupt = length(wk_act{4});
		    ndwnpt = ceil(nupt/2);

		    dwntaxis = [wk_act{4} NaN(1,2*ndwnpt-nupt)];
		    dwntaxis = reshape(dwntaxis,2,ndwnpt);
		    dwntaxis = dwntaxis(1,:);
		    Results.taxis   = dwntaxis;
                    
		    dwnEvol = [wk_Evol{3} NaN(1,2*ndwnpt-nupt)];
		    dwnEvol = reshape(dwnEvol,2,ndwnpt);
		    dwnEvol = dwnEvol(1,:);
		    Results.mass_Evol   = dwnEvol;

		    dwnJtot = [wk_msg{3} NaN(1,2*ndwnpt-nupt)];
		    dwnJtot = reshape(dwnJtot,2,ndwnpt);
		    dwnJtot = dwnJtot(1,:);
		    Results.J_tot = dwnJtot;

		    dwnJdrop = [wk_msg{4} NaN(1,2*ndwnpt-nupt)];
		    dwnJdrop = reshape(dwnJdrop,2,ndwnpt);
		    dwnJdrop = dwnJdrop(1,:);
		    Results.J_drop = dwnJdrop;
		    
                    Results.cG0 = wk_msg{5};
                    Results.Ca0 = wk_msg{6};
                
                    Sbl_results{k}{wkcols(l(idm))} = Results;
                    flag_colpos{k}(wkcols(l(idm))) = true;
                
                    if (pos < length(clock))&&(l(idm)/n_trials*100 >= clock(pos+1))
                        save(['jobid_batch' num2str(i)],...
                            'flag_colpos','Sbl_results','-append');
                        disp(['Job batch ' 'file ' num2str(i) ': ' num2str(k)...
                            '/' num2str(length(jobid)) ' at ' ...
                           num2str(clock(pos+1)) '%'])
                    while (pos < length(clock))&&...
                           (l(idm)/n_trials*100 >= clock(pos+1))
                        pos = pos + 1;
                    end
                    end
                end 
            end
            % Finished this mat
            flag_done(k) = true;
            save(['jobid_batch' num2str(i)],...
                  'flag_colpos','Sbl_results','flag_done',...
                  '-append');
            disp(['Job batch ' 'file ' num2str(i) ': ' num2str(k)...
                  '/' num2str(length(jobid)) ' at ' ...
                  '100%'])
              
            pos = 1;
        end
    end
end
