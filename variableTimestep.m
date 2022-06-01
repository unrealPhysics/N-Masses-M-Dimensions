function [vel,h] = variableTimestep(deltaV,vel,h,binder)
%VARIABLETIMESTEP Returns the timestep, h, after variation and the velocity matrix vel+deltaV
%   Vary the timestep based on the fraction of deltaV div vel according to
%   the value specified in the binder.fractionalAccuracyUpperBound
    fraction = norm(deltaV)/norm(vel);
    if fraction>binder.fractionalAccuracyUpperBound
        h=h/2;
    elseif fraction<binder.fractionalAccuracyUpperBound/2
        h=2*h;
    end
    vel=vel + deltaV;
end

