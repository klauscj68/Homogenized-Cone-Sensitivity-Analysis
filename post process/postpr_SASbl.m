%Post-process the evaluated functionals according to Saltelli02
%  If at some point the number of functionals is changed, you will need to
%  adjust the 'Compute the functionals for base parameter set' and the 
%  'Compute the functional values for all the parameter variations' in
%  addition to changes already made in postpr_functional.
%  'Compute Sensitivity Indices' section has a manual .05 factor that needs
%  to match the job_SApd section Setup parfeval var factor

warning('off','MATLAB:table:ModifiedAndSavedVarnames');
%% Load functionals and initialize
disp('SASbl_Evals.mat and trial_master.mat should be in this folder');
load('trial_master.mat','id_free','n_frparam');
load('Sbl_Evals.mat');
m1 = real(m1);
m2 = real(m2);
for i=1:n_frparam
    ai{i} = real(ai{i}); %#ok<*SAGROW>
    ani{i} = real(ani{i});
end

n_functionals = size(m1,1);
n_trials = size(m1,2);

flags = true(size(m1));
%% Parameter Names
param_names = ...
    {'Rb','Rt','H','cosgamma0','thetain','thetafin','eps0','nu','sig',...
          'flagch1','flagch2','flagch3',...
          'nsez','taglia','tolR','tolangle',...
          'methodcyto1','methodcyto2','methodcyto3','theta','alpha','tolfix','normainf','tfin',...
          'nstept','downsample',...
          'tolstat',...
          'NAv','F',...
          'BcG','BCa',...
          'nuRG','DRst','kR','Rsig',...
          'kGE','DGst','Gsig',...
          'DEst','kE','PDEsig','Betad','Km','kcat',...
          'Esig','Evol','khyd','kcatDIVKm','kst',...
          'utent','kku',...
          'vtent','kkv',...
          'jcGmax','mcG','KcG','fCa',...
          'jexsat','Kex',...
          'amax','amin','mcyc','Kcyc',...
          'nRst0','rate',...
          'flagrestart','rstindex'...
    };
%% Process the MC Means and Variances
% Means
for i=1:n_functionals
    disp(['Functional ' num2str(i)]);
    
    % Very first call also stores ntrials considered not outliers for use
    %  in confidence intervals. It's the same ntrials for all quantities
    %  since threshold is defined in mc_plot
    if i < 12
        [F,G,H,K,...
        mu,sigma,arms,...
        ntrials] = mc_plot(m1(i,flags(i,:)),100); %#ok<*ASGLU>
    else
        [F,G,H,K,...
        mu,sigma,arms,...
        ntrials_chi] = mc_plot(m1(i,flags(i,:)),100,0); %#ok<*ASGLU>
    end
    E.mcest1(i,1) = mu;
    E.mcvar1(i,1) = sigma;
    E.mcconv1(i,1) = arms;
    %pause
    %close all
    
    if i < 12
        [F,G,H,K,...
        mu,sigma,arms] = mc_plot(m2(i,flags(i,:)),100);
    else
        [F,G,H,K,...
        mu,sigma,arms] = mc_plot(m2(i,flags(i,:)),100,0);
    end
    E.mcest2(i,1) = mu;
    E.mcvar2(i,1) = sigma;
    E.mcconv2(i,1) = arms;
    %pause
    %close all
    
end

% Normalize the functionals to have mean 0
shift = .5*(E.mcest1 + E.mcest2);
m1 = m1-shift;
m2 = m2-shift;
for i=1:n_frparam
    ai{i} = ai{i} - shift;
    ani{i} = ani{i} - shift;
end

% Variance
V.sbl = m1.*(m1-m2);
for i=1:n_functionals
    disp(['Functional ' num2str(i)]);
    
    if i < 12
        [F,G,H,K,...
        mu,sigma,arms] = mc_plot(V.sbl(i,flags(i,:)),100);
    else
        [F,G,H,K,...
        mu,sigma,arms] = mc_plot(V.sbl(i,flags(i,:)),100,0);
    end
    
    V.mcest(i,1) = mu;
    V.mcvar(i,1) = sigma;
    V.mcconv(i,1) = arms;
    %pause 
    %close all
end
V.sbl = [];

%% Process Di Indices
for i=n_frparam:-1:1
    Di.sbl1 = m1.*(ai{i} - m2);
    Di.sbl2 = m2.*(ani{i} - m1);
    
    for j=1:n_functionals
        disp(['Di: Param ' param_names{id_free(i)} ' at functional ' num2str(j)]);
        
        if j < 12
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Di.sbl1(j,flags(j,:)),100);
        else
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Di.sbl1(j,flags(j,:)),100,0);
        end
        
        Di.mcest1(j,i) = mu;
        Di.mcvar1(j,i) = sigma;
        Di.mcconv1(j,i) = arms;
     
        %pause
        %close all
    
        if j < 12
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Di.sbl2(j,flags(j,:)),100);
        else
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Di.sbl2(j,flags(j,:)),100,0);
        end
        
        Di.mcest2(j,i) = mu;
        Di.mcvar2(j,i) = sigma;
        Di.mcconv2(j,i) = arms;
        
        %pause
        %close all
    
    end
    
end
Di.sbl1 = [];
Di.sbl2 = [];

%% Process Dtoti Indices
for i=n_frparam:-1:1
    Dtoti.sbl1 = .5*(m1 - ani{i}).^2;
    Dtoti.sbl2 = .5*(m2 - ai{i}).^2;
    
    for j=1:n_functionals
        disp(['Dtoti: Param ' param_names{id_free(i)} ...
              ' at functional ' num2str(j)]);
        
        if j < 12
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Dtoti.sbl1(j,flags(j,:)),100);
        else
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Dtoti.sbl1(j,flags(j,:)),100,0);
        end
        
        Dtoti.mcest1(j,i) = mu;
        Dtoti.mcvar1(j,i) = sigma;
        Dtoti.mcconv1(j,i) = arms;
        %pause
        %close all
    
        if j < 12
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Dtoti.sbl2(j,flags(j,:)),100);
        else
            [F,G,H,K,...
            mu,sigma,arms] = mc_plot(Dtoti.sbl2(j,flags(j,:)),100,0);
        end
        
        Dtoti.mcest2(j,i) = mu;
        Dtoti.mcvar2(j,i) = sigma;
        Dtoti.mcconv2(j,i) = arms;
        %pause
        %close all
    
    end
    
end
Dtoti.sbl1 = [];
Dtoti.sbl2 = [];

%% Process Dtotij Indices
for i=n_frparam:-1:1
    for j=i-1:-1:1
        
        Dtotij.sbl1 = .5*(ai{i} - ai{j}).^2;
        Dtotij.sbl2 = .5*(ani{i} - ani{j}).^2;
    
        for k=1:n_functionals
            disp(['Dtotij: Param ' param_names{id_free(i)} ' ' ...
                  param_names{id_free(j)} ...
                  ' at functional ' num2str(k)]);
        
            if k < 12
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotij.sbl1(k,flags(k,:)),100);
            else
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotij.sbl1(k,flags(k,:)),100,0);
            end
            
            Dtotij.mcest1(k,i,j) = mu;
            Dtotij.mcvar1(k,i,j) = sigma;
            Dtotij.mcconv1(k,i,j) = arms;
            %pause
            %close all
    
            if k < 12
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotij.sbl2(k,flags(k,:)),100);
            else
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotij.sbl2(k,flags(k,:)),100,0);
            end
            
            Dtotij.mcest2(k,i,j) = mu;
            Dtotij.mcvar2(k,i,j) = sigma;
            Dtotij.mcconv2(k,i,j) = arms;
            %pause
            %close all
        end
    
    end
    
end
Dtotij.sbl1 = [];
Dtotij.sbl2 = [];

%% Process Dtotnij Indices
for i=n_frparam:-1:1
    for j=i-1:-1:1
        
        Dtotnij.sbl1 = .5*(ai{j} - ani{i}).^2;
        Dtotnij.sbl2 = .5*(ai{i} - ani{j}).^2;
    
        for k=1:n_functionals
            disp(['Dtotnij: Param ' param_names{id_free(i)} ' ' ...
                  param_names{id_free(j)} ...
                  ' at functional ' num2str(k)]);
        
            if k < 12
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotnij.sbl1(k,flags(k,:)),100);
            else
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotnij.sbl1(k,flags(k,:)),100,0);
            end
            
            Dtotnij.mcest1(k,i,j) = mu;
            Dtotnij.mcvar1(k,i,j) = sigma;
            Dtotnij.mcconv1(k,i,j) = arms;
            %pause
            %close all
    
            if k < 12
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotnij.sbl2(k,flags(k,:)),100);
            else
                [F,G,H,K,...
                mu,sigma,arms] = mc_plot(Dtotnij.sbl2(k,flags(k,:)),100,0);
            end
            
            Dtotnij.mcest2(k,i,j) = mu;
            Dtotnij.mcvar2(k,i,j) = sigma;
            Dtotnij.mcconv2(k,i,j) = arms;
            %pause
            %close all
        end
    
    end
    
end
Dtotnij.sbl1 = [];
Dtotnij.sbl2 = [];

%% Declare functional names for table
ROWNAMES = {'Iact','Ipeak','Irec',...
            'Estact','Estpeak','Estrec','L2',...
            'Tpeak','J+width','Jdark','Jover',...
            'chiIpeak','chiL2','chiTpeak','chiJdark',...
            'chiJover'};

%% Si Sobol Indices
% Write to csv file
fid = fopen('Si.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = Di.mcest1./V.mcest; %Which Saltelli Estimate
ledger = []; %#ok<*NASGU>
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% Si Upper Range Sobol Indices
% Write to csv file
fid = fopen('Si_upper.csv','w');

% Header
header = ['functional'];
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = (Di.mcest1+2*Di.mcvar1/sqrt(ntrials))./(V.mcest-2*V.mcvar/sqrt(ntrials)); %Which Saltelli Estimate
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% Si Lower Range Sobol Indices
% Write to csv file
fid = fopen('Si_lower.csv','w');

% Header
header = ['functional'];
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = (Di.mcest1-2*Di.mcvar1/sqrt(ntrials))./(V.mcest+2*V.mcvar/sqrt(ntrials)); %Which Saltelli Estimate
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% Stoti Sobol Indices
% Write to csv file
fid = fopen('Stoti.csv','w');

% Header
header = ['functional'];
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = Dtoti.mcest1./V.mcest; %Which Saltelli Estimate
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% Stoti Upper Range Sobol Indices
% Write to csv file
fid = fopen('Stoti_upper.csv','w');

% Header
header = ['functional'];
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = (Dtoti.mcest1+2*Dtoti.mcvar1/sqrt(ntrials))./... %Which Saltelli Estimate
       (V.mcest-2*V.mcvar/sqrt(ntrials));
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% Stoti Lower Range Sobol Indices
% Write to csv file
fid = fopen('Stoti_lower.csv','w');

% Header
header = ['functional'];
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = (Dtoti.mcest1-2*Dtoti.mcvar1/sqrt(ntrials))./... %Which Saltelli Estimate
       (V.mcest+2*V.mcvar/sqrt(ntrials));
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% MC V Conv Statistics
% Write to csv file
fid = fopen('Vconv.csv','w');

% Header
header = ['functional,arms'];
header = [header '\n'];
fprintf(fid,header);

% Store average root mean square error
ram = V.mcconv;
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i} ',%5.4f'];
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i));
end

% Close the file
fclose(fid);

%% MC Di Conv Statistics
% Write to csv file
fid = fopen('Diconv.csv','w');

% Header
header = ['functional'];
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store average root mean square error
ram = Di.mcconv1; %Which Saltelli Estimate
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%5.4f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% MC Dtoti Conv Statistics
% Write to csv file
fid = fopen('Dtoticonv.csv','w');

% Header
header = ['functional'];
for i=1:n_frparam
    header = [header ',' param_names{id_free(i)}];
end
header = [header '\n'];
fprintf(fid,header);

% Store average root mean square error
ram = Dtoti.mcconv1; %Which Saltelli Estimate
ledger = [];
for i=1:length(ROWNAMES)
    ledger = [ROWNAMES{i}];
    for j=1:n_frparam
        ledger = [ledger ',%5.4f']; %#ok<*AGROW>
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,ram(i,:));
end

% Close the file
fclose(fid);

%% Stotij Sobol Indices
% Write to csv file
fid = fopen('Stotij.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = Dtotij.mcest1./V.mcest; %Which Saltelli Estimate
ledger = []; %#ok<*NASGU>
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = []; %#ok<*NASGU>
    for i=1:n_frparam
        for j=1:i-1
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
        wram = [wram ram(k,i,j)]; 
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);

%% Stotij Upper Sobol Indices
% Write to csv file
fid = fopen('Stotij_upper.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = (Dtotij.mcest1 + 2*Dtotij.mcvar1/sqrt(ntrials))./... %Which Saltelli Estimate
       (V.mcest - 2*V.mcvar/sqrt(ntrials));
ledger = []; %#ok<*NASGU>
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = [];
    for i=1:n_frparam
        for j=1:i-1
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
        wram = [wram ram(k,i,j)];
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);

%% Stotij Lower Sobol Indices
% Write to csv file
fid = fopen('Stotij_lower.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = (Dtotij.mcest1 - 2*Dtotij.mcvar1/sqrt(ntrials))./... %Which Saltelli Estimate
       (V.mcest + 2*V.mcvar/sqrt(ntrials));
ledger = []; %#ok<*NASGU>
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = [];
    for i=1:n_frparam
        for j=1:i-1
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
        wram = [wram ram(k,i,j)];
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);

%% Sij Sobol Indices
% Write to csv file
fid = fopen('Sij.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = 1-Dtotnij.mcest1./V.mcest; %Which Saltelli Estimate
ledger = []; %#ok<*NASGU>
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = [];
    for i=1:n_frparam
        for j=1:i-1
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
        wram = [wram ram(k,i,j)];
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);

%% Sij Lower Sobol Indices
% Write to csv file
fid = fopen('Sij_lower.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = 1-(Dtotnij.mcest1 + 2*Dtotnij.mcvar1/sqrt(ntrials))./... %Which Saltelli Estimate
       (V.mcest - 2*V.mcvar/sqrt(ntrials));
ledger = []; %#ok<*NASGU>
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = [];
    for i=1:n_frparam
        for j=1:i-1
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
        wram = [wram ram(k,i,j)];
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);

%% Sij Upper Sobol Indices
% Write to csv file
fid = fopen('Sij_upper.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store Sobol indices
ram = 1-(Dtotnij.mcest1 - 2*Dtotnij.mcvar1/sqrt(ntrials))./... %Which Saltelli Estimate
       (V.mcest + 2*V.mcvar/sqrt(ntrials));
ledger = []; %#ok<*NASGU>
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = [];
    for i=1:n_frparam
        for j=1:i-1
        ledger = [ledger ',%4.3f']; %#ok<*AGROW>
        wram = [wram ram(k,i,j)];
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);

%% MC Dtotij Conv Statistics
% Write to csv file
fid = fopen('Dtotijconv.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store average root mean square error
ram = Dtotij.mcconv1; %Which Saltelli Estimate
ledger = [];
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = [];
    for i=1:n_frparam    
        for j=1:i-1
            ledger = [ledger ',%5.4f']; %#ok<*AGROW>
            wram = [wram ram(k,i,j)];
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);

%% MC Dtotnij Conv Statistics
% Write to csv file
fid = fopen('Dtotnijconv.csv','w');

% Header
header = ['functional']; %#ok<*NBRAK>
for i=1:n_frparam
    for j=1:i-1
        % We list like (j,i) because when do ram(k,:,:) matlab cycles along
        % the first coordinate then the second so mimicing that the first
        % coordinate should cycle before the second one
        header = [header ',' param_names{id_free(j)} '/' ...
                             param_names{id_free(i)}];
    end
end
header = [header '\n'];
fprintf(fid,header);

% Store average root mean square error
ram = Dtotnij.mcconv1; %Which Saltelli Estimate
ledger = [];
for k=1:length(ROWNAMES)
    ledger = [ROWNAMES{k}];
    wram = [];
    for i=1:n_frparam    
        for j=1:i-1
            ledger = [ledger ',%5.4f']; %#ok<*AGROW>
            wram = [wram ram(k,i,j)];
        end
    end
    ledger = [ledger '\n'];
    fprintf(fid,ledger,wram);
end

% Close the file
fclose(fid);
%% Common parameters for Box and Whisker Plots
% Define number of (pairs of) parameters to plot
n_boxp = 8;

% Define Caps for Box and Whisker Plots
upcap = 1.25;
lwcap = -.25;
%% Box Whisker Plot for Si
% Extract the names of parameters whose sensitivities are measured
fid = fopen('Si_upper.csv','r');
varnames = textscan(fid,'%s');
varnames = varnames{1}{1};
varnames = split(varnames,',');
fclose(fid);

% Load the ranges
Tu = readtable('Si_upper.csv');
Tl = readtable('Si_lower.csv');

% Loop over functionals to create plots
for i=1:n_functionals
    % Reset from last run
    clear box names
    
    % Setup the box intervals and param names
    box(2,:) = Tu{i,2:end};
    box(1,:) = Tl{i,2:end};
    
    names = varnames(2:end)';
    % Find the greatest quantities
    boxavg = mean(box);
    [~,flag_box] = maxk(boxavg,n_boxp);
    
    %  Cap them before extracting
    flag_cap = (box <= lwcap);
    box(flag_cap) = lwcap;
    flag_cap = (box >= upcap);
    box(flag_cap) = upcap;
    
    box = box(:,flag_box);
    names = names(flag_box);
    
    % Make the plot
    F = myboxplot(box,names);
    
    title([ROWNAMES{i} ' Si indices: 90% Confidence']);
    ylabel('Confidence Intervals');
    
    % Format figure
    F.CurrentAxes.YAxis.Limits = [lwcap upcap];
    F.PaperPosition = [0 0 3.75/F.PaperPosition(3)*F.PaperPosition(3:4)];
    F.PaperSize = F.PaperPosition(3:4);
    F.CurrentAxes.FontName = 'Arial';
    F.CurrentAxes.FontSize = 12;
    
    saveas(F,[ROWNAMES{i} '_Si.fig']);
end
close all

%% Box Whisker Plot for Stoti
% Extract the names of parameters whose sensitivities are measured
fid = fopen('Stoti_upper.csv','r');
varnames = textscan(fid,'%s');
varnames = varnames{1}{1};
varnames = split(varnames,',');
fclose(fid);

% Load the ranges
Tu = readtable('Stoti_upper.csv');
Tl = readtable('Stoti_lower.csv');

% Loop over functionals to create plots
for i=1:n_functionals
    % Reset from last run
    clear box names
    
    % Setup the box intervals and param names
    box(2,:) = Tu{i,2:end};
    box(1,:) = Tl{i,2:end};
    
    names = varnames(2:end)';
    % Find the greatest quantities
    boxavg = mean(box);
    [~,flag_box] = maxk(boxavg,n_boxp);
    
    % Cap them
    flag_cap = (box <= lwcap);
    box(flag_cap) = lwcap;
    flag_cap = (box >= upcap);
    box(flag_cap) = upcap;
    
    box = box(:,flag_box);
    names = names(flag_box);
    
    % Make the plot
    F = myboxplot(box,names);
    
    title([ROWNAMES{i} ' Stoti indices: 90% Confidence']);
    ylabel('Confidence Intervals');
    
    % Format figure
    F.CurrentAxes.YAxis.Limits = [lwcap upcap];
    F.PaperPosition = [0 0 3.75/F.PaperPosition(3)*F.PaperPosition(3:4)];
    F.PaperSize = F.PaperPosition(3:4);
    F.CurrentAxes.FontName = 'Arial';
    F.CurrentAxes.FontSize = 12;
    
    saveas(F,[ROWNAMES{i} '_Stoti.fig']);
end
close all

%% Box Whisker Plot for Stotij
% Extract the names of parameters whose sensitivities are measured
fid = fopen('Stotij_upper.csv','r');
varnames = textscan(fid,'%s');
varnames = varnames{1}{1};
varnames = split(varnames,',');
fclose(fid);

% Load the ranges
Tu = readtable('Stotij_upper.csv');
Tl = readtable('Stotij_lower.csv');

% Loop over functionals to create plots
for i=1:n_functionals
    % Reset from last run
    clear box names
    
    % Setup the box intervals and param names
    box(2,:) = Tu{i,2:end};
    box(1,:) = Tl{i,2:end};
    
    names = varnames(2:end)';
    % Find the greatest quantities
    boxavg = mean(box);
    [~,flag_box] = maxk(boxavg,n_boxp);
    
    %  Cap them before extracting
    flag_cap = (box <= lwcap);
    box(flag_cap) = lwcap;
    flag_cap = (box >= upcap);
    box(flag_cap) = upcap;
    
    box = box(:,flag_box);
    names = names(flag_box);
    
    % Make the plot
    F = myboxplot(box,names);
    
    title([ROWNAMES{i} ' Stotij indices: 90% Confidence']);
    ylabel('Confidence Intervals');
    
    % Format figure
    F.CurrentAxes.YAxis.Limits = [lwcap upcap];
    F.PaperPosition = [0 0 3.75/F.PaperPosition(3)*F.PaperPosition(3:4)];
    F.PaperSize = F.PaperPosition(3:4);
    F.CurrentAxes.FontName = 'Arial';
    F.CurrentAxes.FontSize = 12;
    
    saveas(F,[ROWNAMES{i} '_Stotij.fig']);
end
close all

%% Box Whisker Plot for Sij
% Extract the names of parameters whose sensitivities are measured
fid = fopen('Sij_upper.csv','r');
varnames = textscan(fid,'%s');
varnames = varnames{1}{1};
varnames = split(varnames,',');
fclose(fid);

% Load the ranges
Tu = readtable('Sij_upper.csv');
Tl = readtable('Sij_lower.csv');

% Loop over functionals to create plots
for i=1:n_functionals
    % Reset from last run
    clear box names
    
    % Setup the box intervals and param names
    box(2,:) = Tu{i,2:end};
    box(1,:) = Tl{i,2:end};
    
    names = varnames(2:end)';
    % Find the least quantities
    boxavg = mean(box);
    [~,flag_box] = maxk(boxavg,n_boxp);
    
    % Cap them
    flag_cap = (box <= lwcap);
    box(flag_cap) = lwcap;
    flag_cap = (box >= upcap);
    box(flag_cap) = upcap;
    
    box = box(:,flag_box);
    names = names(flag_box);
    
    % Make the plot
    F = myboxplot(box,names);
    
    title([ROWNAMES{i} ' Sij indices: 90% Confidence']);
    ylabel('Confidence Intervals');
    
    % Format figure
    F.CurrentAxes.YAxis.Limits = [lwcap upcap];
    F.PaperPosition = [0 0 3.75/F.PaperPosition(3)*F.PaperPosition(3:4)];
    F.PaperSize = F.PaperPosition(3:4);
    F.CurrentAxes.FontName = 'Arial';
    F.CurrentAxes.FontSize = 12;
    
    
    saveas(F,[ROWNAMES{i} '_Sij.fig']);
end
close all

%% Save the workspace
%close all
clear F G H

save('Sbl_Analysis.mat');