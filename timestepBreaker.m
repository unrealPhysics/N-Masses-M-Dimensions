function [flag] = timestepBreaker(h,flag)
%TIMESTEPBREAKER Checks whether the timestep has left its bounds and if needed raises the breakFlag

    if log(h) > 10*log(10) || log(h) < -10*log(10)
        disp(h)
        if ~flag.breakFlag
            flag.breakFlag=true;
        end
    end
end

