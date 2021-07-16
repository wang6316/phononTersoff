function [VA] = VA(r,S,De,beta,Re)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%attractive potential VA
VA=-S*De/(S-1)*exp(-beta*sqrt(2/S)*(r-Re));
end

