function A = sbl_dep(M0,M,pointer,param)
%Update any parameters depending on the free parameters
%   M0 is the sample matrix most rows are same as
%   M contains the row that we are substituting into
%   pointer is the index vector of data_set
%   param is the datamat index saying which parameter has changed
%   A is the output resampled matrix

%% Initialize
ntrials = size(M0,2);
A = M0;
A(param,:) = M(param,:);

%% Update by case
if (param == 1)||(param == 2)||(param == 3)
    % cosgamma0 dependended on R_b,R_t,H
    for i=1:ntrials
        indprm = A(1:3,i);
        A(4,i) = fcosgamma0(indprm(1),indprm(2),indprm(3));
    end
end

if param == pointer(8)+3
    % E_sigma dependend on PDE_sigma
    for i=1:ntrials
        indprm = A(pointer(8)+3,i);
        A(pointer(9)+1,i) = fE_sigma(indprm);
    end
end

if (param == pointer(8) + 3)||(param == 8)||(param == 7)
    % E_vol depended on PDE_sigma, nu, epsilon_0
    for i=1:ntrials
        indprm = A([pointer(8)+3,8,7],i);
        A(pointer(9) + 2,i) = fE_vol(indprm(1),indprm(2),indprm(3));
    end
end

if (param == 8) || (param == 7) || (param == pointer(8)+4) ||...
   (param == pointer(8)+3)
   % k_hyd depended on nu,epsilon_0,Beta_dark,PDE_sigma
   for i=1:ntrials
       indprm = A([8,7,pointer(8)+4,pointer(8)+3],i);
       A(pointer(9)+3,i) = fk_hyd(indprm(1),indprm(2),indprm(3),indprm(4));
   end
end

if (param == pointer(9) + 4) || (param == pointer(5) + 1)
    % k_st dependend on kcat_DIV_Km, B_cG
    for i=1:ntrials
        indprm = A([pointer(9) + 4,pointer(5) + 1],i);
        A(pointer(9)+5,i) = fk_st(indprm(1),indprm(2));
    end
end

if (param == pointer(5) + 2)||(param == pointer(12) + 1)
    % j_cG_max depended on B_Ca through its normalizations
    for i=1:ntrials
        oldprm = M0([pointer(5) + 2,pointer(12) + 1,pointer(4) + 2],i);        
        newprm = M([pointer(5) + 2,pointer(12) + 1,pointer(4) + 2],i);
        if param == pointer(12) + 1 %jcgmax case
            jcg_new = newprm(2)*newprm(1)*newprm(3); % j*B_Ca*F
            A(pointer(12)+1,i) = jcg_new/oldprm(1)/oldprm(3);
        else %Bca case
            jcg_old = oldprm(2)*oldprm(1)*oldprm(3);
            indprm = A(pointer(5)+2,i);
            A(pointer(12)+1,i) = jcg_old/indprm/oldprm(3);
        end   
    end
end

if (param == pointer(5) + 2)||(param == pointer(13) + 1) 
    % j_ex_sat depended on B_Ca through its normalizations
    for i=1:ntrials
        oldprm = M0([pointer(5) + 2,pointer(13) + 1,pointer(4) + 2],i);
        newprm = M([pointer(5) + 2,pointer(13) + 1,pointer(4) + 2],i);
        if param == pointer(13) + 1 %jexsat case
            jex_new = newprm(2)*newprm(1)*newprm(3); % j*B_Ca*F
            A(pointer(13)+1,i) = jex_new/oldprm(1)/oldprm(3);
        else %Bca case
            jex_old = oldprm(2)*oldprm(1)*oldprm(3);
            indprm = A(pointer(5)+2,i);
            A(pointer(13)+1,i) = jex_old/indprm/oldprm(3);
        end
    end
end
end

