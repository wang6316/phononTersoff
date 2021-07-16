function [VR] = VR(r,De,S,beta,Re)
%repulsive potential between i,j
%   Detailed explanation goes here
VR=De/(S-1)*exp(-beta*sqrt(2*S)*(r-Re));
end

