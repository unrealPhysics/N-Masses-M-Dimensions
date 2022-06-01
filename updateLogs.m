function [frame,loc,temporal,timeLocked,K,P,T] = updateLogs(frame,mass,pos,vel,loc,temporal,time,timeLocked,K,P,T,flag,constant)
%UPDATELOGS Updates the saved values according to the raised flags and increments the current frame

    frame=frame+1;
    loc(frame,:,:)=pos(:,:); %#ok<*AGROW>
    if flag.temporalInclusionFlag
        temporal(frame)=time;
    end
    if flag.energyFlag
        [K(frame),P(frame),T(frame)] = energyCalcs(mass,pos,vel,constant);
        timeLocked(frame)=time;
    end
end

