function [forcing] = msg_forcing(c_u,c_v,...
                                 E_star,...
                                 M,Msl,Mfo,n_pts,...
                                 flag_ch,...
                                 Sigma_sl,V_cn,Sigma_fo,...
                                 nu,epsilon_0,...
                                 k_hyd,k_st,...
                                 alpha_max,alpha_min,...
                                 m_cyc,K_cyc,...
                                 E_vol,...
                                 B_Ca,F,...
                                 j_cG_max,m_cG,K_cG,f_Ca,...
                                 j_ex_sat,K_ex)                                                                                     
%Compute right hand side of ODE for (cG,Ca)
%   c_u,c_v are the nodal coefficients of resp (cG,Ca) in column form
%   E_star is the nodal coefficients of E_star at present time
%   M,Msl,Mfo are the \int phi_i*phi_j matrices over their resp domains
%   n_pts is the number of nodes in the mesh
%   flag_ch is a logical 3 x 1 saying if channels are at the sliver, discs,
%       or lateral folds
%   Sigma_sl, V_cn, Sigma_fo are respectively the surface areas of the
%       sliver, volume of the cone interior and surface area of folds. Can 
%        be obtained by double summing the matching M matrix.
%   nu is the ratio of interdiscal space to discs of increment epsilon_0
%   k_hyd, k_st are the hydrolysis rates of cG by dark and activated E
%   alpha_max, alpha_min,m_cyc,K_cyc are biochemical parameters of cyclase
%       synthesis of cG
%   E_vol is the total volumic concentration of all E in the cone
%   B_Ca is the buffering capacity of calcium in the cytoplasm
%   F is Faraday's constant
%   j_cG_max,m_cG,K_cG,f_Ca are biochemical parameters for cG-gated current
%   j_ex_sat,K_ex are biochemical parameters of exchanger current

%% cG hydrolysis
% Compute forcing contribution of dark and light cG hydrolysis before \int
F11 = (nu/(1+nu))*(k_hyd/2-k_st)*(E_star.*c_u);


%% cG synthesis by cyclase
% Compute the forcing contribution of cG synthesis by cyclase before \int
F12 = (nu/(1+nu))*(...
                   alpha_min + ...
                   (alpha_max - alpha_min)./(1 + (c_v./K_cyc).^m_cyc)...
                  );

%% Integrate and aggregate the F1 term for cG right hand side
F1 = M*(F11+F12);

%% Initialize channel normalizations
eta0 = .5*nu*epsilon_0;
Sigma_Ch = flag_ch(1)*Sigma_sl + ...
           flag_ch(2)*(nu/(1+nu))*V_cn/eta0 + ...
           flag_ch(3)*Sigma_fo;
       
%% Preliminary computation of current densities 
% cG-gated current
J_cG = j_cG_max/Sigma_Ch*(c_u.^m_cG)./(K_cG.^m_cG + c_u.^m_cG);
J_cG = (.5*f_Ca)*J_cG;

% Exchanger current
J_ex = j_ex_sat/Sigma_Ch*(c_v./(K_ex+c_v));
J_ex = J_ex;

%% Ca2+ influx at sliver
if flag_ch(1) == true    
    F21 = Msl*(J_cG - J_ex);
else
    F21 = zeros(n_pts,1);
end

%% Ca2+ influx at discs
if flag_ch(2) == true
    F22 = (nu/(1+nu))*(1/eta0)*(M*(J_cG - J_ex));
else
    F22 = zeros(n_pts,1);
end

%% Ca2+ influx at folds
if flag_ch(3) == true
    F23 = (nu/(1+nu))*(Mfo*(J_cG - J_ex));
else
    F23 = zeros(n_pts,1);
end

%% Aggregate F2 term for Ca2+ right hand side
F2 = F21 + F22 + F23;

%% Create forcing term
forcing = [F1;...
           F2];

end


