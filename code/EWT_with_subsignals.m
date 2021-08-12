function [vital_signal,subsignals_br,subsignals_hr]=EWT_with_subsignals(vital_signal,band,bounds)
    vital_signal=normalization(vital_signal,band);   
    f=vital_signal;
    f=f';
    t=0:length(f)-1;
    params.SamplingRate = 20; %put -1 if you don't know the sampling rate
    % Choose the wanted global trend removal (none,plaw,poly,morpho,tophat,opening)
    params.globtrend = 'none';
    params.degree=3; % degree for the polynomial interpolation
    % Choose the wanted regularization (none,gaussian,average,closing)
    params.reg = 'none';
    params.lengthFilter = 10;
    params.sigmaFilter = 1.5;
    % Choose the wanted detection method (locmax,locmaxmin,
    % adaptive,adaptivereg,scalespace)
    params.detect = 'vitalsigns';
    params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans
    params.N = 5; % maximum number of bands
    params.completion = 0; % choose if you want to force to have params.N modes
                           % in case the algorithm found less ones (0 or 1)
    params.init_bounds = bounds*length(vital_signal)/20;
    params.log=0;

    %% We perform the empirical transform and its inverse
    [ewt,mfb,boundaries]=EWT1D(f,params);
    subsignals_br=[];
    subsignals_hr=[];
    if size(bounds,2)==4
        for j=1:size(ewt,1)
            tmpewt=ewt{j};
            tmp=abs(fftshift(fft(tmpewt,2^16)));
            [~,maxI]=max(tmp(32768:end));
            freq=maxI/32768*10*60;
            if bounds(1)*60<freq && freq<bounds(2)*60
                subsignals_br=[subsignals_br,tmpewt];
            end

            if bounds(3)*60<freq && freq<bounds(4)*60
                subsignals_hr=[subsignals_hr,tmpewt];
            end  
        end
    else       
        for j=1:size(ewt,1)
            tmpewt=ewt{j};
            tmp=abs(fftshift(fft(tmpewt,2^16)));
            [~,maxI]=max(tmp(32768:end));
            freq=maxI/32768*10*60;
            if bounds(1)*60<freq && freq< bounds(2)*60
                subsignals_br=[subsignals_br,tmpewt];
            end          
        end
        subsignals_hr=[subsignals_hr,ewt{end-1}];
    end
end

function sig = normalization(raw_sig,band)
%     butter_sig=smoothdata(raw_sig,'gaussian',3);
    max1=max(raw_sig);
    min1=min(raw_sig);
    sig=(raw_sig-min1)/(max1-min1)-0.5;    
end