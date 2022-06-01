function [loc,timeLocked,K,P,T] = initialiseSmalls(pos)
%INITIALISESMALLS Initialises a bunch of values in a single line

loc(1,:,:)=pos(:,:);
timeLocked(1)=0;
K(1)=0;
P(1)=0;
T(1)=0;
end

