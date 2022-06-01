H=newplot;
reset(H);
clear;
format compact
H=newplot;
reset(H);

N=3;
M=2;
pos(N,M)=0;
vel(N,M)=0;
acc(N,M)=0;
accN(N,M)=0;
mass=ones(N,1);
G=1;

h=0.01;
steps=20000;
skipVal=round(steps/100);
firstLoop=true;

pos(1,:)=[1,0];
vel(1,:)=[0,1];
pos(2,:)=[-1,0];
vel(2,:)=[0,-1];
pos(3,:)=[0,0];
vel(3,:)=[0,0];
mass=[1;1;2];

for i=1:steps
    loc(i,:,:)=pos(:,:);

    if rem(i,skipVal)==0
        axis on;
        plot(loc(:,1,1),loc(:,1,2),'g-',loc(:,2,1),loc(:,2,2),'b-',...
             loc(:,3,1),loc(:,3,2),'r-')
        axis equal;
        title(i)
        drawnow
    end

    for body1=1:N
        for body2=(body1+1):N
            r((body2-body1)+(body1-1)*(N-body1+1))=norm(pos(body1,:)-pos(body2,:));
        end
    end

    if firstLoop % fix and do things
        for body1=1:N
            for body2=(body1+1):N
                r((body2-body1)+(body1-1)*(N-body1+1))=norm(pos(body1,:)-pos(body2,:));
            end
        end
        acc=zeros(N,M);
        for body1=1:N
            for body2=(body1+1):N
                accVec=G*mass(body2)/r(body1,body2)^3 * (pos(body1,:)-pos(body2,:));
                acc(body1,:)=acc(body1,:)-accVec;
                acc(body2,:)=acc(body2,:)+accVec;
            end
        end
        firstLoop=false;
    end

    pos=pos + h*vel + 0.5*h*h*acc;

    for body1=1:N
        for body2=(body1+1):N
            r(body1,body2)=norm(pos(body1,:)-pos(body2,:));
        end
    end
    
    accN=zeros(N,M);
    for body1=1:N
        for body2=(body1+1):N
            accVec=G*mass(body2)/r(body1,body2)^3 * (pos(body1,:)-pos(body2,:));
            accN(body1,:)=accN(body1,:)-accVec;
            accN(body2,:)=accN(body2,:)+accVec;
        end
    end

    vel=vel + 0.5*h*(acc+accN);

    acc=accN;

end



