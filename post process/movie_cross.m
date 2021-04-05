function [MV] = movie_cross(u,taxis,twin,...
                            z,...
                            var_name)
%Make a movie of given solution
% u is a n_pts x n_taxis sol array we wish to plot
% taxis is a row vector of the time points columns of u were sampled at
% twin is a 1 x 2 indicating the starting and finishing times for the movie
%  stands for t window
% z is either empty, a 1 x 1, or 1 x 2.  Respectively, these tell the code
%  to average over the whole cone, at the single level z, or the sections
%  spanned in z by z(1),z(2)
% varname is the string used to title the plot

%% Parameters
% Number of movie frames. Default playback is 30 frames/s
nframes = 360;

% Setup colorbar limits
ti = twin(1);
tf = twin(2);
j_ti = find(taxis <= ti,1,'last');%  greatest lower bound
j_tf = find(taxis >= tf,1);%  least upper bound

CLimit = [min(min(u(:,j_ti:j_tf))) max(max(u(:,j_ti:j_tf)))];

%% Times taking movie frames
dt = tf - ti;
tframes = ti:dt/(nframes-1):tf;
if size(tframes,2) ~= nframes
    tframes = [tframes tf];
    
    assert(size(tframes,2) == nframes,...
           'Did not achieve correct number of frames');
end

%% Create Movie
% Preallocate structure which stores movie frames
MV(nframes) = struct('cdata',[],'colormap',[]);

% Loop over and save frames
WBar = waitbar(0,'Creating Movie Frames');

for i=1:nframes
    % Get current time point
    t = tframes(i);
    
    % Build frame, format, and save to movie
    plot_cross(u,taxis,t,...
               z,...
               var_name);
           
    F = gcf;
    F.CurrentAxes.CLim = CLimit;
    
    MV(i) = getframe(F);
    
    % Close figure so can make next figure
    close(F);
    
    % Update waitbar
    waitbar(i/nframes,WBar);
end
close(WBar);

%% Save movie to mpg4 format
WBar = waitbar(0,'Saving Movie to MPG4');

V = VideoWriter(var_name,'MPEG-4');
open(V);

for i=1:nframes
    writeVideo(V,MV(i));
    waitbar(i/nframes,WBar);
end
close(WBar)

close(V);

end

