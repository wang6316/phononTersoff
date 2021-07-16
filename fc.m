function [fc] = fc(r,R,D)
%cutoff function base on equation(8) 
%   Detailed explanation goes here
if (r<=(R-D))
    fc=1;
elseif(r<(R+D))
    fc=0.5-0.5*sin(pi*(r-R)/(2*D));
else
    fc=0;
end

end

