function [K,P,T] = energyCalcs(mass,pos,vel,constant)
%ENERGYCALCS Calculates the kinetic, potential and total energies of the system

    P=0;
    K=0;
    for body1=1:constant.N
        for body2=(body1+1):constant.N
            P=P - constant.G * mass(body1)*mass(body2)/(norm(pos(body1,:)-pos(body2,:)));
        end
        K=K+0.5*mass(body1)*norm(vel(body1,:))^2;
    end

    T=P+K;
end

