function [vital_signal]=EWT_on_row(raw_data)
    pure_data=raw_data(:,1:end-3)-mean(raw_data(:,1:end-3));
    recn=pure_data;
    for k=1:size(pure_data,1)
        close all
        f=pure_data(k,:);
        f=f';
        t=0:length(f)-1;
        params.SamplingRate = -1; %put -1 if you don't know the sampling rate
        % Choose the wanted global trend removal (none,plaw,poly,morpho,tophat,opening)
        params.globtrend = 'poly';
        params.degree=3; % degree for the polynomial interpolation
        % Choose the wanted regularization (none,gaussian,average,closing)
        params.reg = 'none';
        params.lengthFilter = 10;
        params.sigmaFilter = 1.5;
        % Choose the wanted detection method (locmax,locmaxmin,
        % adaptive,adaptivereg,scalespace)
        params.detect = 'locmaxminf';
        params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans
        params.N = 5; % maximum number of bands
        params.completion = 0; % choose if you want to force to have params.N modes
%         % Perform the detection on the log spectrum instead the spectrum
        params.log=0;
%         InitBounds = params.InitBounds;

        %% We perform the empirical transform and its inverse
        % compute the EWT (and get the corresponding filter bank and list of 
        % boundaries)
        [ewt,mfb,boundaries]=EWT1D(f,params);
        ens=[];
        for num=1:size(ewt,1)
            ens=[ens sum(ewt{num}.^2)];
        end
        [~,maxI]=max(ens);
        recn(k,:)=ewt{maxI};
    end
    
    maxindexs=[];
    maxAms=[];
    for i=1:size(recn)
        [maxAm,maxindex]=max(recn(i,:));
        maxindexs=[maxindexs,maxindex];
        maxAms=[maxAms,recn(i,maxindex)];
    end
    vital_signal=maxAms;    
end