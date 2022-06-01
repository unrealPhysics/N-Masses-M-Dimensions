function [acc] = fullAccGen(effectiveMass,pos,constant)
%FULLACCGEN Generates the full acceleration matrix

    influenceMat=zeros(constant.N,constant.N,constant.M);

    for body1=1:(constant.N-1)
        for body2=(body1+1):constant.N
            vector=(pos(body2,:)-pos(body1,:))/(norm(pos(body1,:)-pos(body2,:))^3);
            influenceMat(body1,body2,:)=effectiveMass(body2)*vector;
            influenceMat(body2,body1,:)=-effectiveMass(body1)*vector;
        end
    end
    acc=reshape(sum(influenceMat,2),constant.N,constant.M);
end

