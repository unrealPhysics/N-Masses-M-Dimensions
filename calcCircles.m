% sqrt(20*(2*(1/(41)^2*(40/41))))
% 
% (20/(sqrt(20*(2*(1/(41)^2*(40/41))))))/(9/(sqrt(9*(1/(18)^2+2/(41)^2*(9/41)))))

rSmall=1;
% rTotal=sqrt(4*rBiiig^2+rSmall^2);
% periodBig=(rBiiig/(sqrt(rBiiig*(2*(1/rTotal^2)*(rBig/rTotal)))));
% periodSmall=(rSmall/(sqrt(rSmall*(1/(4*rSmall^2)+(2/rTotal^2)*(rSmall/rTotal)))));
% periodBig=10*periodSmall;
% 
% (rBiiig/(sqrt(rBiiig*(2*(1/rTotal^2)*(rBig/rTotal)))))=10*(rSmall/(sqrt(rSmall*(1/(4*rSmall^2)+(2/rTotal^2)*(rSmall/rTotal)))));

% for n=7:1000
%     rTotal=2*(n)^(2/3);
%     % rTotal=(4*(100*pi^2-2))^(1/3)*rSmall;
%     % rBiiig=sqrt((4*(800*pi^2-2))^(2/3)-1)*rSmall/2;
%     % rTotal=sqrt(4*rBiiig^2+rSmall^2);
%     rBiiig=sqrt(rTotal^2-rSmall^2)/2;
%     vBiiig=sqrt(rBiiig*(2*rTotal^(-2)*2*rBiiig/rTotal));
%     vSmall=sqrt(rSmall*(1/(4*rSmall^2)+(2/rTotal^2)*(rSmall/rTotal)));
%     periodBiiig=(rBiiig/(sqrt(rBiiig*(2*(1/rTotal^2)*(2*rBiiig/rTotal)))));
%     periodSmall=(rSmall/(sqrt(rSmall*(1/(4*rSmall^2)+(2/rTotal^2)*(rSmall/rTotal)))));
%     periodFrac=periodBiiig/periodSmall;
%     if rem(periodFrac,1)==0
%         disp(periodFrac)
%         disp(periodFrac)
%     end
% end

% (rBiiig/(sqrt(rBiiig*(2*(1/rTotal^2)*(2*rBiiig/rTotal)))))/(rSmall/(sqrt(rSmall*(1/(4*rSmall^2)+(2/rTotal^2)*(rSmall/rTotal)))))

% v=sqrt(r*a)
% vBiiig=sqrt(rBiiig*(2*rTotal^(-2)*2*rBiiig/rTotal))


periodsPerPeriod=2
rTotal=rSmall*(4*((2*periodsPerPeriod)^2-2))^(1/3);
rBiiig=sqrt(rTotal^2-rSmall^2)/2;
vBiiig=sqrt(rBiiig*(2*rTotal^(-2)*2*rBiiig/rTotal));
vSmall=sqrt(rSmall*(1/(4*rSmall^2)+(2/rTotal^2)*(rSmall/rTotal)));
periodBiiig=(rBiiig/(sqrt(rBiiig*(2*(1/rTotal^2)*(2*rBiiig/rTotal)))));
periodSmall=(rSmall/(sqrt(rSmall*((1/(4*rSmall^2))+(2/rTotal^2)*(rSmall/rTotal)))));
periodFrac=periodBiiig/periodSmall;
% 
% periodBiiig=sqrt(rTotal^3)/2;
% 
% periodSmall=sqrt(rSmall*rTotal^3/((rTotal^3/(4*rSmall^2))+(2*rSmall)));
% % periodSmall=rSmall/((1/(4*rSmall^2))+(2*rSmall/rTotal^3));
% periodFrac2=0.5*sqrt((rTotal^3/(4*rSmall^3)) + 2);
% rTotal=rSmall*(4*((2*periodFrac2)^2-2))^(1/3);
% %     rTotal=rSmall*(4*(2*periodFrac2)^2-2)^(1/3)


