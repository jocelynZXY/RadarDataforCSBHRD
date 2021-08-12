function [br,hr] = corr_vital_signs(subsignals_br1,subsignals_br2, subsignals_hr1, subsignals_hr2)
    coffs=zeros(size(subsignals_br1,2),size(subsignals_br2,2));
    for k=1:size(subsignals_br1,2) 
        for j=1:size(subsignals_br2,2)
              tmpcorr=xcorr(subsignals_br1(:,k),subsignals_br2(:,j));
              coffs(k)=max(tmpcorr);
        end
    end
    [x,y]=find(coffs==max(max(coffs)));
    [crossf,~]=xcorr(subsignals_br1(:,x),subsignals_br2(:,y));

    tmp=abs(fftshift(fft(crossf,2^16)));
    [~,maxI]=max(tmp(32768:end));
    br=maxI/32768*10*60;

    coffs=zeros(size(subsignals_hr1,2),1);
    
    for k=1:size(subsignals_hr1,2)
        tmp=corrcoef(subsignals_hr1(:,k),subsignals_hr2);
        coffs(k)=tmp(1,2);
    end
    
    [~,maxI]=find(coffs==max(coffs));
     en_hr=[];
    for num=1:size(subsignals_hr1,2)
        if num~=maxI
            en_hr=[en_hr sum(subsignals_hr1(:,num).^2)];
        else
            en_hr=[en_hr 0];
        end
    end
    
    [~,maxI]=max(en_hr);
    tmp=abs(fftshift(fft(subsignals_hr1(:,maxI),2^16)));
    [~,maxI]=max(tmp(32768:end));
    hr=maxI/32768*10*60;
end

