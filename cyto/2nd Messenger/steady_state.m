% soluzione numerica del sistema non lineare che descrive la soluzione stazionaria
function [u_ss,v_ss]=steady_state(R_b,R_t,H,theta_in,theta_fin,epsilon_0,nu,...
                                  flag_ch,...
                                  tol_stat,...
                                  B_Ca,F,...
                                  PDE_sigma,k_hyd,...
                                  u_tent,...
                                  v_tent,...
                                  j_cG_max,m_cG,K_cG,f_Ca,...
                                  j_ex_sat,K_ex,...
                                  alpha_max,alpha_min,m_cyc,K_cyc)

%fprintf('\nSoluzione dello steady-state\n');

%% Run iterative flux balance

% dato iniziale
x0(1)=u_tent;
x0(2)=v_tent;

% risolve
FS = @(x) fun(x, epsilon_0, nu, ...
    k_hyd, PDE_sigma, alpha_max, alpha_min, m_cyc, K_cyc, ...
    B_Ca, F, j_cG_max, f_Ca, m_cG, K_cG, j_ex_sat, K_ex);

x = fsolve(FS, x0, optimset('TolFun', tol_stat,...
                            'Display','off')); 

% mette la soluzione in u_ss e v_ss
u_ss=real(x(1));
v_ss=real(x(2));

%disp(['Steady State: [cGMP] = ' num2str(u_ss) ' uM '...
%                'and [Ca2+] = ' num2str(v_ss) ' uM']);
                                


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun(x, epsilon_0, nu, ...
    k_hyd, PDE_sigma, alpha_max, alpha_min, m_cyc, K_cyc, ...
    B_Ca, F, j_cG_max, f_Ca, m_cG, K_cG, j_ex_sat, K_ex)
% estrae le incognite u,v da x
u=x(1);
v=x(2);

% output
y=zeros(1,2);

% prima equazione:
% bilancio dei flussi del cGMP sulle facce dei dischi
y(1)=k_hyd*PDE_sigma*u-...
          (alpha_min+(alpha_max-alpha_min)*K_cyc^m_cyc/(K_cyc^m_cyc+v^m_cyc))*(1/2*nu*epsilon_0);



% seconda equazione:
% bilancio dei flussi del Ca2+ sulla membrana plasmatica
y(2)=j_ex_sat*v/(K_ex+v)-j_cG_max*f_Ca/(2)*u^m_cG/(K_cG^m_cG+u^m_cG);
