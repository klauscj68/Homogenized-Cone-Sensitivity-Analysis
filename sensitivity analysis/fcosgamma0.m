function [cosgamma0] = fcosgamma0(R_b,R_t,H)
%Compute cosgamma0 for given geometry

h = abs(R_t*H/...
       (R_b-R_t));
   
if isinf(h)
    cosgamma0 = 1;
else
    cosgamma0 = h/sqrt(h^2+R_t^2);
end

end

