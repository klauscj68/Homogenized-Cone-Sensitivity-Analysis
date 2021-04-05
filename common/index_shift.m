function [shift] = index_shift(n_pts,n_sez)
%Index n_sez copies of pts with their copy value.
%   shift is an array like [1 .. 1 2 .. 2 ... n_sez ... n_sez] 
%   there are as many duplicates as there are n_pts

shift = zeros(1,n_pts);
shift(1) = 1;
shift = repmat(shift,1,n_sez);
shift = cumsum(shift);

end

