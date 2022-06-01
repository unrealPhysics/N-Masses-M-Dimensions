function [pos,vel,mass,N,M,G,axis] = initialPresetter(preset)
%INITIALPRESETTER Returns a predefined group of initial values according to the input value, preset

    switch preset
        case 1
            N=3;
            M=3;
            G=1;
            mass=ones(1,N);

            pos(1,:)=[ 1, 0, 0];
            vel(1,:)=[ 0,-1, 0];
            pos(2,:)=[-1, 0, 0];
            vel(2,:)=[ 0, 1, 0];
            pos(3,:)=[ 0, 0, 0];
            vel(3,:)=[ 0, 0, 0];

            axis=[-1,1,-1,1,-1,1];
        case 2
            N=6;
            M=3;
            G=1;
            mass=ones(1,N);

            pos(1,:)=[ 0, 0, 1];
            pos(2,:)=[ 0, 0,-1];
            pos(3,:)=[ 0, 1, 0];
            pos(4,:)=[ 0,-1, 0];
            pos(5,:)=[ 1, 0, 0];
            pos(6,:)=[-1, 0, 0];
            
            vel(1,:)=[ 0, 1, 0];
            vel(2,:)=[ 0,-1, 0];
            vel(3,:)=[ 1, 0, 0];
            vel(4,:)=[-1, 0, 0];
            vel(5,:)=[ 0, 0, 1];
            vel(6,:)=[ 0, 0,-1];

            axis=[-1,1,-1,1,-1,1];
        case 3
            N=5;
            M=2;
            G=1;

            mass=ones(1,N);

            pos(1,:)=[ 1, 0];
            pos(2,:)=[-1, 0];
            pos(3,:)=[ 0, 0];
            pos(4,:)=[ 0, 1];
            pos(5,:)=[ 0,-1];
            
            vel(1,:)=[ 0,-1];
            vel(2,:)=[ 0, 1];
            vel(3,:)=[ 0, 0];
            vel(4,:)=[ 1, 0];
            vel(5,:)=[-1, 0];

            axis=[-1,1,-1,1];
        case 4
            N=5;
            M=3;
            G=1;

            mass=2*ones(1,N);

            pos(1,:)=[ 1, 0,-1];
            pos(2,:)=[-1, 0,-1];
            pos(3,:)=[ 0, 0, 1];
            pos(4,:)=[ 0, 1,-1];
            pos(5,:)=[ 0,-1,-1];
            
            vel(1,:)=[ 0,-2, 0];
            vel(2,:)=[ 0, 2, 0];
            vel(3,:)=[ 0, 0, 0];
            vel(4,:)=[ 2, 0, 0];
            vel(5,:)=[-2, 0, 0];

            mass(3)=2*4;

            axis=2*[-1,1,-1,1,-1,1];
        case 5
            N=3;
            M=3;
            G=1;

            mass=ones(1,N);

            pos(1,:)=[ 1, 0, 0];
            pos(2,:)=[ 0, 1, 0];
            pos(3,:)=[ 0, 0, 0];
            
            vel(1,:)=[ 0, 0, 1];
            vel(2,:)=[ 0, 0,-1];
            vel(3,:)=[ 0, 0, 0];

            axis=[-1,1,-1,1,-1,1];
        case 6
            N=4;
            M=4;
            G=1;

            mass=ones(1,N);

            pos(1,:)=[-1,0,0,0];
            vel(1,:)=[0,0,0,1];
            pos(2,:)=[1,0,0,0];
            vel(2,:)=[0,0,0,-1];
            pos(3,:)=[0,0,0,1];
            vel(3,:)=[1,0,0,0];
            pos(4,:)=[0,0,0,-1];
            vel(4,:)=[-1,0,0,0];

            axis=[-1,1,-1,1,-1,1];
        case 7
            % distances in AU
            % speeds in km/s
            % masses in M_Earth
%             % G currently hecked
            N=9;
            M=2;
            G=1;

            mass=ones(1,N);

            pos(1,:)=[0,0];
            vel(1,:)=[0,0];
            mass(1)=332900; %sun

            pos(2,:)=[0.4,0];
            vel(2,:)=[0,47.36];
            mass(2)=0.055; %mercury

            pos(3,:)=[0.66666666,0];
            vel(3,:)=[0,35.02];
            mass(3)=0.815; %venus

            pos(4,:)=[1,0];
            vel(4,:)=[0,29.78];
            mass(4)=1; %earth

            pos(5,:)=[1.5,0];
            vel(5,:)=[0,24.007];
            mass(5)=0.107; %mars

            pos(6,:)=[5.2,0];
            vel(6,:)=[0,12.07];
            mass(6)=318; %jupiter

            pos(7,:)=[9.5,0];
            vel(7,:)=[0,9.68];
            mass(7)=95; %saturn

            pos(8,:)=[19.2,0];
            vel(8,:)=[0,6.8];
            mass(8)=14; %uranus

            pos(9,:)=[30.1,0];
            vel(9,:)=[0,5.43];
            mass(9)=17; %neptune

            %Conversion to SI units

            G=6.674*10^(-11)*G;
            mass=5.97237*10^24*mass; %M_Earth to kg
            pos=150*10^9*pos; %AU to m
            vel=10^(3)*vel; %km/s to m/s

            axis=norm(pos(9,:))*[-1,1,-1,1];
        case 8
            N=2;
            M=2;
            G=1;

            mass=3*ones(1,N);

            pos(1,:)=[0,1];
            vel(1,:)=[1,0];
            pos(2,:)=[0,-1];
            vel(2,:)=[-1,0];

            axis=[-1,1,-1,1];
        case 9
            N=4;
            M=3;
            G=1;

            mass=3*ones(1,N);
            mass(3)=0;
            mass(4)=0;

            pos(1,:)=[0,1,0];
            vel(1,:)=[1,0,0];
            pos(2,:)=[0,-1,0];
            vel(2,:)=[-1,0,0];
            pos(3,:)=[0,1.1,0];
            vel(3,:)=[1,0,4];
            pos(4,:)=[0,-1.1,0];
            vel(4,:)=[-1,0,-4];

            axis=[-1,1,-1,1,-1,1];
        case 10
            N=3;
            M=4;
            G=1;

            mass=ones(1,N);
            mass(3)=sum(mass(1:2));

            rSmall=1;
            periodsPerPeriod=0.99; % must be greater than 0.75
            rTotal=rSmall*(4*((2*periodsPerPeriod)^2-2))^(1/3);
            rBiiig=sqrt(rTotal^2-rSmall^2)/2;
            vBiiig=sqrt(rBiiig*(2*rTotal^(-2)*2*rBiiig/rTotal));
            vSmall=sqrt(rSmall*(1/(4*rSmall^2)+(2/rTotal^2)*(rSmall/rTotal)));

            pos(1,:)=rSmall*[0,0,0,1];
            vel(1,:)=vSmall*[0,0,1,0];
            pos(2,:)=-pos(1,:);
            vel(2,:)=-vel(1,:);
            pos(3,:)=rBiiig*[1,0,0,0];
            vel(3,:)=vBiiig*[0,1,0,0];
            pos(1,:)=pos(1,:)-pos(3,:);
            pos(2,:)=pos(2,:)-pos(3,:);
            vel(1,:)=vel(1,:)-vel(3,:);
            vel(2,:)=vel(2,:)-vel(3,:);

            axis=(floor(rBiiig)+1)*[-1,1,-1,1,-1,1];

%             2*(sqrt(3^2+4^2)^(-4))*3+(3)^(-3)
%             sqrt(4*(2*(sqrt(3^2+4^2)^(-3))*4))
% 
%             v=sqrt(R*G*(M1/r1^2+M2/r2^2))
%             v=sqrt(R*(M1/r1^2+M2/r2^2))
% 
%             v=sqrt(3*(1/(3)^2+2/(5)^2*3/5))
%             0.6908931418
% %             v=sqrt(4*(2/(5)^2*4/5))
% %             0.5059644256
% 
%             v=sqrt(R_orbitpoint*accel)
%             accel=G*(m1/d1^2+m2/d2^2)
%             asmall=(1/(18)^2+2/(41)^2*(9/41))
%             sqrt(9*(1/(18)^2+2/(41)^2*(9/41)))
%             alorg=(2*(1/(41)^2*(40/41)))
%             sqrt(20*(2*(1/(41)^2*(40/41))))
% %             
        case 11
            N=3;
            M=3;
            G=1;

            mass=2*ones(1,N);
            mass(3)=sum(mass(1:2));

            pos(1,:)=[0,0,0.1];
            vel(1,:)=2.24*[0,1,0]; %should be circles
            pos(2,:)=-pos(1,:);
            vel(2,:)=-vel(1,:);
            pos(3,:)=[1,0,0];
            vel(3,:)=[0,1,0];
            pos(1,:)=pos(1,:)-pos(3,:);
            pos(2,:)=pos(2,:)-pos(3,:);
            vel(1,:)=vel(1,:)-vel(3,:);
            vel(2,:)=vel(2,:)-vel(3,:);

            axis=[-1,1,-1,1,-1,1];
        case 12
            N=2;
            M=2;
            G=1;

            mass=zeros(1,N);
            mass(1)=2;

            pos(1,:)=[0,0];
            vel(1,:)=[1,0];
            pos(2,:)=[0,1];
            vel(2,:)=[1,1];

            axis=[-1,1,-1,1];
        case 13
            N=5;
            M=3;
            G=1;

            mass=0.1*ones(1,N);
            mass(3)=1;

            pos(1,:)=[ 1, 0, 0];
            vel(1,:)=[ 0, 1, 0];
            pos(2,:)=[-1, 0, 0];
            vel(2,:)=[ 0,-1, 0];
            pos(3,:)=[ 0, 0, 0];
            vel(3,:)=[ 0, 0, 0];
            pos(4,:)=[ 0, 1, 0];
            vel(4,:)=[ 0, 0, 1];
            pos(5,:)=[ 0,-1, 0];
            vel(5,:)=[ 0, 0,-1];

            axis=2*[-1,1,-1,1,-1,1];
        otherwise
            disp('Invalid preset!')
            disp('Loading default preset (1)...')
            [pos,vel,mass,N,M,G,axis] = initialPresetter(1);
    end
end

