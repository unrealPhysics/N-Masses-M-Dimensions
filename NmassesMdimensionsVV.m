function NmassesMdimensionsVV()
%NMASSESMDIMENSIONSVV Simulates N-point masses in M-dimensions using a Velocity-Verlet integrator

format compact
figure(2)

% Set choice flags
% Essentially is the settings of the input/output
flag.tripleDimPlotFlag=true;
flag.energyFlag=true;
flag.equalityFlag=true;
flag.temporalInclusionFlag=false;
flag.enforcedBoundsFlag=true;
flag.usePresets=true;
preset=10;

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

acc(constant.N,constant.M)=0;

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
flag.firstLoopFlag=true;
flag.breakFlag=false;
flag.finalPlotFlag=false;

% Initialise the non-program related values
effectiveMass=constant.G*mass;
[loc,timeLocked,K,P,T] = initialiseSmalls(pos);

% Start the timer
tInitial=tic;

for i=1:binder.steps
    if flag.firstLoopFlag
        [acc]=fullAccGen(effectiveMass,pos,constant);
        flag.firstLoopFlag=false;
    end
    
    if rem(i-1,binder.iterFrame)==0
        [frame,loc,temporal,timeLocked,K,P,T] = updateLogs(frame,mass,pos,vel,loc,temporal,time,timeLocked,K,P,T,flag,constant);
    end
    time=time+h;

    % Plot
    currTime = tocPlot(i,frame,tInitial,time,currTime,pos,loc,temporal,timeLocked,K,P,T,constant,flag,binder);

    % Round the coordinates to kill epsilon errors
    pos(:)=round(pos(:),10,'significant');
    vel(:)=round(vel(:),10,'significant');
    acc(:)=round(acc(:),10,'significant');

    % Velocity-Verlet
    pos=pos + h*vel + 0.5*h*h*acc;

    % Generate the next acceleration matrix
    [accN]=fullAccGen(effectiveMass,pos,constant);

    % Calculate deltaV via Velocity-Verlet Method
    deltaV=0.5*h*(acc+accN);

    % Variable timestep & vel=vel+deltaV
    [vel,h] = variableTimestep(deltaV,vel,h,binder);

    % Set the current acceleration matrix equal to the next acceleration
    % matrix
    acc=accN;

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
% kerotage=norm(pos(3,:));
% for i=1:constant.N
%     kelemel=norm(reshape(loc(i,3,:),constant.M,1));
%     if kelemel>kerotage
%         kerotage=kelemel;
%     end
% end
% disp(kerotage)
end
