function [rate, flag_tdep] = flash(pts_eul,t,...
                                   rate,flag_query)
%Compute light flash intensity across the nodes of mesh
%   pts_eul is 3 x n_pts coords of nodes across Eul coord mesh of cone 
%   t is time in seconds
%   rate is the fraction of unactivated R that is activated by the flash
%       per second. As output, it is a column vector matching the indexing 
%       of pts_eul. As input, it is the value specified in data_set. When a
%       number, the flash is uniform and constant in space and time per the
%       input number. When it is NaN, the script defaults the function
%       defined by isnan(rate) == true below.
%   flag_query is an optional argument and if true, the program returns NaN
%       for rate and the value of flag_dep. If flag_query is true, dummy
%       values may be passed for pts_eul and t.  The script will not use
%       them.
%   flag_tdep is true if the flash varies in time and false if it is if it 
%       is a constant.  It needs to be set by the user below, and is 
%       important because this flag affects how the global matrices
%       are assembled. Obviously, it must match the given rule for rate of
%       activation.
%   The user must manually edit lines 24 and section "Space-time Varying
%       Illumination" if giving a custom flash function.

%% User Declare if Space-time varying flash is t dependent
flag_tdep = true;

% Proofread and if rate passed was number, then we are not t dependent
if ~isnan(rate)
    flag_tdep = false;
end

%% Sample Rate Over Mesh
if nargin == 3
    switch isnan(rate)
        case true
            %% Space-time Varying Illumination
            % Initialize variables
            n_pts = size(pts_eul,2);
        
            % Compute rate of activation at time point
            rate = repmat(100*exp(-t),n_pts,1);
                
        case false
            %% Constant Illimunation
            % The rate is constant and uniform in space and time
            n_pts = size(pts_eul,2);
            rate = repmat(rate,n_pts,1);
        
            % Set flag_tdep
            flag_tdep = false;
        
    end
end


