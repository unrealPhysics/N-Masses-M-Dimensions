function NmassesMdimensionsOptimised()
%NMASSESMDIMENSIONS Summary of this function goes here
%   Detailed explanation goes here
format compact
% H=newplot;
% reset(H);
figure(2)

global globals %#ok<*GVMIS>
flags.tripleDimPlotFlag=true;
flags.energyFlag=true;
flags.usePresets=true;
flags.enforcedBoundsFlag=false;
preset=3 %#ok<NOPRT> 

if flags.usePresets
    [pos,vel,mass,globals.N,globals.M,G]=initialPresetter(preset);
else
    globals.N=6;
    globals.M=3;
    G=1;
    pos(globals.N,globals.M)=0;
    vel(globals.N,globals.M)=0;
end

acc(globals.N,globals.M)=0;

h=0.0001;
steps=80000;
skipVal=round(steps/100);
fractionalAccuracyUpperBound=0.002;
time=0;
kel=0;
flags.firstLoopFlag=true;
flags.breakFlag=false;

if globals.M<3 && flags.tripleDimPlotFlag
    flags.tripleDimPlotFlag=false;
    disp('tripleDimPlotFlag has been set to false due to initial dimensions less than three');
end
tInitial=tic;
effectiveMass=G*mass;
for i=1:steps
    if flags.firstLoopFlag
%         massMat=defineMassMat(mass,G);
        [acc]=fullAccGen(effectiveMass,pos);
        flags.firstLoopFlag=false;
    end
    
    loc(i,:,:)=pos(:,:); %#ok<*AGROW>
    time=time+h;
    if flags.energyFlag
        [K(i),P(i),T(i)] = energyCalcs(mass,pos,vel,G);
        timeLocked(i)=time;
    end
    

    if rem(i,skipVal)==0
        axis on;
        if flags.energyFlag
            subplot(1,2,1)
            if flags.tripleDimPlotFlag
                plot3(pos(:,1),pos(:,2),pos(:,3),'o',loc(:,:,1),loc(:,:,2),loc(:,:,3))
            else
                if flags.enforcedBoundsFlag
                    plot(pos(:,1),pos(:,2),'o',loc(:,:,1),loc(:,:,2))
                    axis([-4.65e12,4.65e12,-4.65e12,4.65e12])
                else
                    plot(pos(:,1),pos(:,2),'o',loc(:,:,1),loc(:,:,2))
                end
            end
            axis equal;
            subplot(1,2,2)
            plot(timeLocked,K,'r-',timeLocked,P,'b-',timeLocked,T,'k-')
        else
            if flags.tripleDimPlotFlag
                plot3(pos(:,1),pos(:,2),pos(:,3),'o',loc(:,:,1),loc(:,:,2),loc(:,:,3))
            else
                if flags.enforcedBoundsFlag
                    plot(pos(:,1),pos(:,2),'o',loc(:,:,1),loc(:,:,2))
                    axis([-4.65e12,4.65e12,-4.65e12,4.65e12])
                else
                    plot(pos(:,1),pos(:,2),'o',loc(:,:,1),loc(:,:,2))
                end
            end
            axis equal;
        end
        kel=kel+1;
        runTime(kel)=toc(tInitial);
        title([time i])
        drawnow
        if flags.breakFlag
            break
        end
    end

    pos(:)=round(pos(:),10,'significant');
    vel(:)=round(vel(:),10,'significant');
    acc(:)=round(acc(:),10,'significant');
    pos=pos + h*vel + 0.5*h*h*acc;
%     pos(:)=round(pos(:),10,'significant');

    [accN]=fullAccGen(effectiveMass,pos);

    deltaV=0.5*h*(acc+accN);
    fraction = norm(deltaV)/norm(vel);
    if fraction>fractionalAccuracyUpperBound
        h=h/2;
    elseif fraction<fractionalAccuracyUpperBound/2
        h=2*h;
    end
    vel=vel + deltaV;

    acc=accN;

    if log(h) > 10*log(10) || log(h) < -10*log(10)
        disp(h)
        if flags.breakFlag==false
            flags.breakFlag=true;
        end
    end

end
tFinal=toc(tInitial);
sprintf('Time elapsed %f',tFinal)
figure(3)
hold on
plot(runTime)
hold off
% plot(T)
disp(max(T)-min(T))
end

function [massMat] = defineMassMat(mass,G)
    global globals
    for body1=1:globals.N
        for body2=1:globals.N
            if body1==body2
                massMat(body1,body2)=0;
            else
                massMat(body1,body2)=mass(body2);
            end
        end
    end
    massMat=G*massMat;
end

function [acc] = fullAccGen(mass,pos)
    global globals
    
%     for body1=1:globals.N
%         for body2=1:body1-1
%             masslessInfluenceVecMat(body2,:)=(pos(body2,:)-pos(body1,:))/(norm(pos(body1,:)-pos(body2,:))^3);
%         end
%         for body2=(body1+1):globals.N
%             masslessInfluenceVecMat(body2,:)=(pos(body2,:)-pos(body1,:))/(norm(pos(body1,:)-pos(body2,:))^3);
%         end
%         acc(body1,:)=massMat(body1,:)*masslessInfluenceVecMat(:,:);
%     end

    influenceMat=zeros(globals.N,globals.N,globals.M);

    for body1=1:(globals.N-1)
        for body2=(body1+1):globals.N
            vector=(pos(body2,:)-pos(body1,:))/(norm(pos(body1,:)-pos(body2,:))^3);
            influenceMat(body1,body2,:)=mass(body2)*vector;
            influenceMat(body2,body1,:)=-mass(body1)*vector;
        end
    end
    acc=reshape(sum(influenceMat,2),globals.N,globals.M);
end

function [K,P,T] = energyCalcs(mass,pos,vel,G)
    global globals

    P=0;
    K=0;
    for body1=1:globals.N
        for body2=(body1+1):globals.N
            P=P - G*mass(body1)*mass(body2)/(norm(pos(body1,:)-pos(body2,:)));
        end
        K=K+0.5*mass(body1)*norm(vel(body1,:))^2;
    end

    T=P+K;
end
