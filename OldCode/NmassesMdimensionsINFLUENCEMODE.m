function NmassesMdimensionsINFLUENCEMODE()
%NMASSESMDIMENSIONS Summary of this function goes here
%   Detailed explanation goes here
format compact
H=newplot;
reset(H);

global N %#ok<*GVMIS> 
N=5;
global M
M=3;
shell=1;

% for shell=[1,-1]
% for bodySwapper=[1,N]
G=1;
pos(N,M)=0;
vel(N,M)=0;
acc(N,M)=0;

h=0.0001;
steps=80000;
skipVal=round(steps/100);
frame=1;
time=0;
fps=300;
firstLoop=true;
breakFlag=false;

spacialTranslation=[0,0,0];

% if bodySwapper==N
%     shell=-shell;
    pos(1,:)=[ 1, 0, 0]+spacialTranslation;
    vel(1,:)=[ 0, 1, 0];
    pos(2,:)=[-1, 0, 0]+spacialTranslation;
    vel(2,:)=[ 0,-1, 0];
    pos(3,:)=[ 0, 0, 0]+spacialTranslation;
    vel(3,:)=[ 0, 0, 0];
    pos(4,:)=[ 0, 1, 0]+spacialTranslation;
    vel(4,:)=[ 0, 0, 1];
    pos(5,:)=[ 0,-1, 0]+spacialTranslation;
    vel(5,:)=[ 0, 0,-1];
% else
%     pos(1,:)=shell*[ 1, 0, 0]-spacialTranslation;
%     vel(1,:)=shell*[ 0, 1, 0];
%     pos(2,:)=shell*[-1, 0, 0]-spacialTranslation;
%     vel(2,:)=shell*[ 0,-1, 0];
%     pos(3,:)=shell*[ 0, 0, 0]-spacialTranslation;
%     vel(3,:)=shell*[ 0, 0, 0];
%     pos(4,:)=shell*[ 0, 1, 0]-spacialTranslation;
%     vel(4,:)=shell*[ 0, 0, 1];
%     pos(5,:)=shell*[ 0,-1, 0]-spacialTranslation;
%     vel(5,:)=shell*[ 0, 0, 1];
% end

% pos(indexer(1,bodySwapper),:)=shell*[ 1, 0, 0];
% vel(indexer(1,bodySwapper),:)=shell*[ 0, 1, 1];
% pos(indexer(2,bodySwapper),:)=shell*[-1, 0, 0];
% vel(indexer(2,bodySwapper),:)=shell*[ 0,-1, 1];
% pos(indexer(3,bodySwapper),:)=shell*[ 0, 0, 0];
% vel(indexer(3,bodySwapper),:)=shell*[ 0, 0, 1];
% pos(indexer(4,bodySwapper),:)=shell*[ 0, 1, 0];
% vel(indexer(4,bodySwapper),:)=shell*[-1, 0, 1];
% pos(indexer(5,bodySwapper),:)=shell*[ 0,-1, 0];
% vel(indexer(5,bodySwapper),:)=shell*[ 1, 0, 1];
mass=[0.1;0.1;1;0.1;0.1];
    for i=1:steps
        if firstLoop
            massMat=defineMassMat(mass,G);
            [acc]=fullAccGen(massMat,pos);
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
                            specialBodies(specialCount)=(N+1)*body1+body2;
                        end
                    end
                end
            end
            firstLoop=false;
        end
%         pos(:,1)=pos(:,1)-pos(3,1);
%         pos(:,2)=pos(:,2)-pos(3,2);
        loc(i,:,:)=pos(:,:); %#ok<*AGROW> 
    %     for body=1:N
    %         loggedTimes(i)=time;
    %         EK(i,body)=mass(body)*norm(vel(body,:))^2;
    %         EP(i,body)=-mass(body)*norm(acc(body,:))*norm(pos(body,:)-pos(3,:));
    %         ET(i,body)=EK(i,body)+EP(i,body);
    %     end
        time=time+h;
    
        if time*fps>frame
            frame=frame+1;
    %         loc(frame,:,:)=pos(:,:);
        end
        if rem(i,skipVal)==0
    %         subplot(2,1,1)
            axis on;
            plot3(loc(:,1,1),loc(:,1,2),loc(:,1,3),'g-',loc(:,2,1),loc(:,2,2),loc(:,2,3),'b-',...
                 loc(:,3,1),loc(:,3,2),loc(:,3,3),'r-',...
                 loc(:,4,1),loc(:,4,2),loc(:,4,3),'m-',loc(:,5,1),loc(:,5,2),loc(:,5,3),'c-',...
                 pos(1,1),pos(1,2),pos(1,3),'go',pos(2,1),pos(2,2),pos(2,3),'bo',...
                 pos(3,1),pos(3,2),pos(3,3),'ro',...
                 pos(4,1),pos(4,2),pos(4,3),'mo',pos(5,1),pos(5,2),pos(5,3),'co')
            axis equal;
            title([time i])
            drawnow
    %         subplot(2,1,2)
    %         plot(loggedTimes(:),EK(:,1),'g-',loggedTimes(:),EK(:,4),'b-',...
    %              loggedTimes(:),EK(:,3),'r-',...
    %              loggedTimes(:),EK(:,4),'m-',loggedTimes(:),EK(:,5),'c-',...
    %              loggedTimes(:),EP(:,1),'g.',loggedTimes(:),EP(:,4),'b.',...
    %              loggedTimes(:),EP(:,3),'r.',...
    %              loggedTimes(:),EP(:,4),'m.',loggedTimes(:),EP(:,5),'c.',...
    %              loggedTimes(:),ET(:,1),'go',loggedTimes(:),ET(:,4),'bo',...
    %              loggedTimes(:),ET(:,3),'ro',...
    %              loggedTimes(:),ET(:,4),'mo',loggedTimes(:),ET(:,5),'co')
    %         drawnow
            if breakFlag
                break
            end
        end
    
        pos=pos + h*vel + 0.5*h*h*acc;
    
        [accN]=fullAccGen(massMat,pos);
    
        if specialCount~=0
            [specialCount,specialBodies] = instabilityCheck(accN,specialBodies,specialCount,i);
        end
    
        deltaV=0.5*h*(acc+accN);
        fraction = norm(deltaV)/norm(vel);
        if fraction>0.001
            h=h/2;
        elseif fraction<0.0005
            h=2*h;
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

%     keroSAVE(((3-shell)/2),:)=loc(:);
%     if shell==1
%         keroSAVE(1,:,:,:)=loc(:,:,:);
%     else
%         keroSAVE(2,:,:,:)=loc(:,:,:);
%     end
%     clear loc
% end
% kevonel=reshape((keroSAVE(1,:)+keroSAVE(2,:))/2,i,N,M);
% plot(kevonel(:,1,1),kevonel(:,1,2),'g-',kevonel(:,2,1),kevonel(:,2,2),'b-',...
%                  kevonel(:,3,1),kevonel(:,3,2),'r-',...
%                  kevonel(:,4,1),kevonel(:,4,2),'m-',kevonel(:,5,1),kevonel(:,5,2),'c-')
% for temel=0:(N-1)
%     kevonel(:,temel+1,:)=reshape((keroSAVE(1,:,deindexer(1,temel),:)+keroSAVE(2,:,deindexer(2,temel),:)+keroSAVE(3,:,deindexer(3,temel),:)+keroSAVE(4,:,deindexer(4,temel),:)+keroSAVE(5,:,deindexer(5,temel),:))/N,i,1,M);
% end
% 
% kevonel(:,1,:)=reshape((keroSAVE(1,:,1,:)+keroSAVE(2,:,1,:))/2,i,1,M);
% kevonel(:,2,:)=reshape((keroSAVE(1,:,2,:)+keroSAVE(2,:,2,:))/2,i,1,M);
% kevonel(:,3,:)=reshape((keroSAVE(1,:,3,:)+keroSAVE(2,:,3,:))/2,i,1,M);
% kevonel(:,4,:)=reshape((keroSAVE(1,:,4,:)+keroSAVE(2,:,4,:))/2,i,1,M);
% kevonel(:,5,:)=reshape((keroSAVE(1,:,5,:)+keroSAVE(2,:,5,:))/2,i,1,M);
% 
% 
% 
% plot3(kevonel(:,1,1),kevonel(:,1,2),kevonel(:,1,3),'g-',kevonel(:,2,1),kevonel(:,2,2),kevonel(:,2,3),'b-',...
%                  kevonel(:,3,1),kevonel(:,3,2),kevonel(:,3,3),'r-',...
%                  kevonel(:,4,1),kevonel(:,4,2),kevonel(:,4,3),'m-',kevonel(:,5,1),kevonel(:,5,2),kevonel(:,5,3),'c-')

%kevonel(:,3,:) defines a locked body, thus apply pseudolock???

end

function [massMat] = defineMassMat(mass,G)
    global N
    for body1=1:N
        for body2=1:N
            if body1==body2
                massMat(body1,body2)=0;
            else
                massMat(body1,body2)=mass(body2);
            end
        end
    end
    massMat=G*massMat;
end

function [acc] = fullAccGen(massMat,pos)
    global N
    
    for body1=1:N
        for body2=1:body1-1
            masslessInfluenceVecMat(body2,:)=(pos(body2,:)-pos(body1,:))/(norm(pos(body1,:)-pos(body2,:))^3);
        end
        for body2=(body1+1):N
            masslessInfluenceVecMat(body2,:)=(pos(body2,:)-pos(body1,:))/(norm(pos(body1,:)-pos(body2,:))^3);
        end
        acc(body1,:)=massMat(body1,:)*masslessInfluenceVecMat(:,:);
    end
    acc(3,:)=0*acc(3,:);
end

function [specialCountNew,specialBodiesNew] = instabilityCheck(acc,specialBodies,specialCount,iter)
    global N

    specialCountNew=specialCount;
    for i=1:specialCount
        if specialBodies(i)==0
            %do nothing
        elseif specialBodies(i) <= N %be fixed body
            zero=acc(specialBodies(i),:);
            if norm(zero)~=0
                fprintf('\nWARNING! WARNING! INSTABILITY DETECTED!\nBODY %d IS HAS BECOME UNFIXED!\nITERATION BE %d\n',specialBodies(i),iter)
                specialCountNew=specialCountNew-1;
                specialBodies(i)=0;
                disp(zero)
            end
        else
            body1=floor(specialBodies(i)/(N+1));
            body2=specialBodies(i)-(N+1)*body1;
            zero=acc(body1,:)+acc(body2,:);
            if norm(zero)~=0
                fprintf('\nWARNING! WARNING! INSTABILITY DETECTED!\nBODIES %d and %d HAVE BECOME UNLINKED!\nITERATION BE %d\n',body1,body2,iter)
                specialCountNew=specialCountNew-1;
                specialBodies(i)=0;
                disp(zero)
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

function index = indexer(normalVal,swapVal)
    global N
    sum=normalVal+swapVal;
    index=rem(sum,N);
    if index==0
        index=N;
    end
end

function deindex = deindexer(index,swapVal)
    global N
    deindex=index-swapVal;
    while deindex > N
        deindex=deindex-N;
    end
    while deindex < 1
        deindex=deindex+N;
    end
end


