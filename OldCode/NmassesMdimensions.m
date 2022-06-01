function NmassesMdimensions()
%NMASSESMDIMENSIONS Summary of this function goes here
%   Detailed explanation goes here
format compact
H=newplot;
reset(H);

global N %#ok<*GVMIS> 
N=5;
global M
M=3;
global G
G=1;
global bit
bit=true;
pos(N,M)=0;
vel(N,M)=0;
acc(N,M)=0;

h=0.0005;
steps=10000;
skipVal=round(steps/100);
frame=1;
time=0;
fps=300;
firstLoop=true;
breakFlag=false;

pos(1,:)=[ 1, 0, 0];
vel(1,:)=[ 0, 1, 1];
pos(2,:)=[-1, 0, 0];
vel(2,:)=[ 0,-1, 1];
pos(3,:)=[ 0, 0, 0];
vel(3,:)=[ 0, 0, 1];
pos(4,:)=[ 0, 1, 0];
vel(4,:)=[-1, 0, 1];
pos(5,:)=[ 0,-1, 0];
vel(5,:)=[ 1, 0, 1];
mass=[1;1;1;1;1];

for i=1:steps
    loc(i,:,:)=pos(:,:); %#ok<*AGROW> 
    time=time+h;

    if time*fps>frame
        frame=frame+1;
%         loc(frame,:,:)=pos(:,:);
    end
    if rem(i,skipVal)==0
        axis on;
        plot(loc(:,1,1),loc(:,1,2),'g-',loc(:,2,1),loc(:,2,2),'b-',...
             loc(:,3,1),loc(:,3,2),'r-',...
             loc(:,4,1),loc(:,4,2),'m-',loc(:,5,1),loc(:,5,2),'c-',...
             pos(1,1),pos(1,2),'go',pos(2,1),pos(2,2),'bo',...
             pos(3,1),pos(3,2),'ro',...
             pos(4,1),pos(4,2),'mo',pos(5,1),pos(5,2),'co')
        axis equal;
        title([time i])
        drawnow
        if breakFlag
            break
        end
    end

    if firstLoop
        r=defineR(pos);
        [acc]=accelGen(mass,r,pos);
        specialCount=0;
        for body1=1:N
            if norm(acc(body1,:))==0
                specialCount=specialCount+1;
                fprintf('Body %d is fixed!\n',body1)
                specialBodies(specialCount)=body1;
            else
                for body2=(body1+1):N
                    if norm(acc(body1,:)+acc(body2,:))==0
                        specialCount=specialCount+1;
                        fprintf('Body %d and body %d are linked!\n',body1,body2)
                        specialBodies(specialCount)=N+index(body1,body2);
                    end
                end
            end
        end
        firstLoop=false;
    end

    pos=pos + h*vel + 0.5*h*h*acc;

    r=defineR(pos);

%     if i==18
%         fprintf('aaaa')
%     end
    
    [accN]=accelGen(mass,r,pos);

    if specialCount~=0
        [specialCount,specialBodies,accN] = instabilityCheck(accN,specialBodies,specialCount,i);
    end

    deltaV=0.5*h*(acc+accN);
    fraction = norm(deltaV)/norm(vel);
    if fraction>0.002
        h=h/2;
    elseif fraction<0.001
%         h=2*h;
    end
    vel=vel + deltaV;

    acc=accN;

    if log(h) > 10*log(10) || log(h) < -10*log(10)
        disp(h)
        if breakFlag==false
            breakFlag=true;
        end
    end

end
end

function index = index(body1,body2)
    global N
    index = body2 + N*body1 - 0.5*(body1^2+body1+2*N);
end

function [r] = defineR(pos)
    global N
    for body1=1:(N-1)
        for body2=(body1+1):N
            r(index(body1,body2))=norm(pos(body1,:)-pos(body2,:));
        end
    end
end

function [acc] = accelGen(mass,r,pos)
    global G
    global N
    global M
    acc=zeros(N,M);
    for body1=1:(N-1)
        for body2=(body1+1):N
            accVec1=G*mass(body2)/r(index(body1,body2))^3 * (pos(body1,:)-pos(body2,:));
            accVec2=G*mass(body1)/r(index(body1,body2))^3 * (pos(body2,:)-pos(body1,:));
%             if(accVec1*mass(body1)~=accVec2*mass(body2))
%                 fprintf('AAAAAAAAAAAAAAAAAAA')
%             end
            acc(body1,:)=acc(body1,:)-accVec1;
            acc(body2,:)=acc(body2,:)-accVec2;
        end
    end
end

function [specialCountNew,specialBodiesNew,acc] = instabilityCheck(acc,specialBodies,specialCount,iter)
    global N
    global bit

    specialCountNew=specialCount;
    for i=1:specialCount
        if specialBodies(i)==0
            %do nothing
        elseif specialBodies(i) <= N %be fixed body
            zero=acc(specialBodies(i),:);
            if norm(zero)~=0
                disp(zero)
                fprintf('\nWARNING! WARNING! INSTABILITY DETECTED!\nBODY %d IS HAS BECOME UNFIXED!\nITERATION BE %d\n',specialBodies(i),iter)
                specialCountNew=specialCountNew-1;
                specialBodies(i)=0;
            end
        else
            [body1,body2] = deIndex(specialBodies(i)-N);
            zero=acc(body1,:)+acc(body2,:);
            if norm(zero)~=0
                disp(zero)
                fprintf('\nWARNING! WARNING! INSTABILITY DETECTED!\nBODIES %d and %d HAVE BECOME UNLINKED!\nITERATION BE %d\n',body1,body2,iter)
%                 fprintf('\n\nAGGRESSIVE RESTABILISATION MODE ACTIVE!\n\n')
%                 if bit
%                     acc(body1,:)=-acc(body2,:);
%                 else 
%                     acc(body2,:)=-acc(body1,:);
%                 end
                    specialCountNew=specialCountNew-1;
                    specialBodies(i)=0;
            end
        end
    end
    if specialCount~=specialCountNew
        j=0;
        specialBodiesNew(1)=0;
        for i=1:specialCount
            if specialBodies(i)~=0
                j=j+1;
                specialBodiesNew(j)=specialBodies(i);
            end
        end
    else
        specialBodiesNew=specialBodies;
    end
end

function [body1,body2] = deIndex(index)
    global N

    trueBody1=0;
    for body1=1:N
        if (body1+1 + N*body1 - 0.5*(body1^2+body1+2*N))>index && trueBody1==0 
            trueBody1=body1-1;
        end
    end
    body1=trueBody1;
    body2=index-(N*body1 - 0.5*(body1^2+body1+2*N));
end
