function NmassesMdimensionsRK()
%NMASSESMDIMENSIONSRK Simulates N-point masses in M-dimensions using a fourth order Runge-Kutta integrator

format compact
figure(1)

% Set choice flags
% Essentially is the settings of the input/output
flag.tripleDimPlotFlag=false;
flag.energyFlag=true;
flag.equalityFlag=true;
flag.temporalInclusionFlag=false;
flag.enforcedBoundsFlag=false;
flag.usePresets=true;
preset=8;

% Enable the use of preset values
if flag.usePresets
    fprintf('Running Preset %d\n',preset)
    [pos,vel,mass,constant.N,constant.M,constant.G,binder.axis]=initialPresetter(preset);
else
    % Default Simulation
    constant.N=3;
    constant.M=3;
    constant.G=1;
    mass=ones(1,N);
    
    pos(1,:)=[ 1, 0, 0];
    vel(1,:)=[ 0,-1, 0];
    pos(2,:)=[-1, 0, 0];
    vel(2,:)=[ 0, 1, 0];
    pos(3,:)=[ 0, 0, 0];
    vel(3,:)=[ 0, 0, 0];

    binder.axis=[-1,1,-1,1,-1,1];
end

if constant.M<3 && flag.tripleDimPlotFlag
    flag.temporalInclusionFlag=true;
    disp('temporalInclusionFlag has been set to true due to tripleDimPlotFlag being set to true while M<3')
end

% Set the required values
% Essentially is the settings of the program
h=0.0001;
binder.steps=400000;
binder.iterFrame=10;
binder.fps=30;
binder.fractionalAccuracyUpperBound=0.002;

% Initialise the program related values
time=0;
frame=0;
currTime=0;
temporal=0;
flag.breakFlag=false;
flag.finalPlotFlag=false;

% Initialise the non-program related values
effectiveMass=constant.G*mass;
[loc,timeLocked,K,P,T] = initialiseSmalls(pos);

% Start the timer
tInitial=tic;

for i=1:binder.steps
    if rem(i-1,binder.iterFrame)==0
        [frame,loc,temporal,timeLocked,K,P,T] = updateLogs(frame,mass,pos,vel,loc,temporal,time,timeLocked,K,P,T,flag,constant);
    end
    time=time+h;
    
    % Plot
    currTime = tocPlot(i,frame,tInitial,time,currTime,pos,loc,temporal,timeLocked,K,P,T,constant,flag,binder);

    % Round the coordinates to kill epsilon errors
    pos(:)=round(pos(:),10,'significant');
    vel(:)=round(vel(:),10,'significant');

    % Fourth order Runge-Kutta
    [pos,deltaV] = RK4(effectiveMass,pos,vel,h,constant);

    % Variable timestep & vel=vel+deltaV
    [vel,h] = variableTimestep(deltaV,vel,h,binder);

    % Check whether the timestep is out of bounds and if so break the
    % program
    [flag] = timestepBreaker(h,flag);

    if flag.breakFlag
        break
    end

end
flag.finalPlotFlag=true;
[frame,loc,temporal,timeLocked,K,P,T] = updateLogs(frame,mass,pos,vel,loc,temporal,time,timeLocked,K,P,T,flag,constant);
tocPlot(i,frame,tInitial,time,currTime,pos,loc,temporal,timeLocked,K,P,T,constant,flag,binder);
tFinal=toc(tInitial);
fprintf('Time elapsed %f seconds\n',tFinal)
end

function [pos,deltaV] = RK4(mass,pos,vel,h,constant)
    % fullAccGen functions as a derivative function
    F1acc=fullAccGen(mass,pos,constant);
    F1vel=vel;
    F2acc=fullAccGen(mass,pos+0.5*h*F1vel,constant);
    F2vel=vel+0.5*h*F1acc;
    F3acc=fullAccGen(mass,pos+0.5*h*F2vel,constant);
    F3vel=vel+0.5*h*F2acc;
    F4acc=fullAccGen(mass,pos+h*F3vel,constant);
    F4vel=vel+h*F3acc;
    pos=pos+h/6*(F1vel+2*F2vel+2*F3vel+F4vel);
    deltaV=h/6*(F1acc+2*F2acc+2*F3acc+F4acc);
end