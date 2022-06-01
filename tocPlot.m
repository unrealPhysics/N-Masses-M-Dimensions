function currTime = tocPlot(i,frame,tInitial,time,currTime,pos,loc,temporal,timeLocked,K,P,T,constant,flag,binder)
%TOCPLOT Plots values as defined by the raised flags and updates the current time

    clockTime=toc(tInitial);
    if clockTime>currTime/binder.fps || flag.finalPlotFlag
        currTime=floor(clockTime*binder.fps)+2;
        axis on;
        if flag.equalityFlag
            axis equal
        end
        if flag.energyFlag
            subplot(1,2,1)
            if flag.tripleDimPlotFlag
                if flag.temporalInclusionFlag
                    plot3(pos(:,1),pos(:,2),temporal(frame)*ones(constant.N,1),'o',loc(:,:,1),loc(:,:,2),temporal.*ones(constant.N,1))
                else
                    plot3(pos(:,1),pos(:,2),pos(:,3),'o',loc(:,:,1),loc(:,:,2),loc(:,:,3))
                    if flag.enforcedBoundsFlag
                        axis(binder.axis)
                    end
                end
            else
                plot(pos(:,1),pos(:,2),'o',loc(:,:,1),loc(:,:,2))
                if flag.enforcedBoundsFlag && length(binder.axis)==4
                    axis(binder.axis)
                end
            end
            if flag.equalityFlag
                axis equal
            end
            subplot(1,2,2)
            plot(timeLocked,K,'r-',timeLocked,P,'b-',timeLocked,T,'k-')
        else
            if flag.tripleDimPlotFlag
                if flag.temporalInclusionFlag
                    plot3(pos(:,1),pos(:,2),temporal(frame)*ones(constant.N,1),'o',loc(:,:,1),loc(:,:,2),temporal.*ones(constant.N,1))
                else
                    plot3(pos(:,1),pos(:,2),pos(:,3),'o',loc(:,:,1),loc(:,:,2),loc(:,:,3))
                    if flag.enforcedBoundsFlag
                        axis(binder.axis)
                    end
                end
            else
                plot(pos(:,1),pos(:,2),'o',loc(:,:,1),loc(:,:,2))
                if flag.enforcedBoundsFlag && length(binder.axis)==4
                    axis(binder.axis)
                end
            end
            if flag.equalityFlag
                axis equal
            end
        end
        title([time i])
        drawnow
    end
end

