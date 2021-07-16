function [Vij] = Vij(r,R,D,bij,De,S,beta,Re)
%VIJ Summary of this function goes here
%   Detailed explanation goes here
Vij=fc(r,R,D)*(VR(r,De,S,beta,Re)+bij*VA(r,S,De,beta,Re));

end

