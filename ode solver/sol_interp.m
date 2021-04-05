function [sol_val] = sol_interp(sol,tpts,t)
%Compute the value of sol at time t
%   sol is n_pts x n_tpts array of sol spline to be interpolated 
%   tpts are the time points columns of sol are sampled at
%   t is the time at which we want to interpolate the sol

%% Find interval containing t
index = find(t <= tpts,1);
assert(~isempty(index),'t is not in interval given by tpts');
index = index - 1;

%% Compute Eval
% Interpolate value between the ends of interval
if index == 0
    % If was at very first time instant, then copy that sol value
    sol_val = sol(:,1);
else
    dt = tpts(index+1)-tpts(index);
    sol_val = (tpts(index+1)-t)/dt*sol(:,index) + ...
              (t-tpts(index))/dt*sol(:,index+1);
end

end

