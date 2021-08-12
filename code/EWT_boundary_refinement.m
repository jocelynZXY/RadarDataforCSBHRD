function bound = EWT_boundary_refinement(init_bounds,f)
    sub_br=boundary_refinement(f,[init_bounds(1),init_bounds(2)],3);
    sub_hr=boundary_refinement(f,[init_bounds(3),init_bounds(4)],5);
    bound=[sub_br,sub_hr];
end

function bound=boundary_refinement(f,band,N)
    f2=f;
    locmin=max(f2)*ones(size(f2));
    locmax=zeros(size(f2));

    minIndex=[];
    for i=2:length(f)-1
        if ((f(i-1)<f(i)) && (f(i)>f(i+1)))
            locmax(i)=f(i);
        end
        if ((f2(i-1)>f2(i)) && (f2(i)<f2(i+1)))
             locmin(i)=f2(i);
             minIndex=[minIndex,i];
        end
    end
    
    for i=1:size(band,2)
        d=abs(minIndex-band(i));
        [~,tmpI]=min(d);
        band(i)=minIndex(tmpI);
    end    
    tmpf=f(band(1):band(2));
    locmax=locmax(band(1):band(2));
    locmin=locmin(band(1):band(2));

    f=tmpf;
    f2=f;   
    
    indexs=find(locmax>0);
    validN=length(indexs);
    if validN+1<N
        N=validN+1;
    end
       
    if N~=-1
        N=N-1;
        [lmax,Imax]=sort(locmax,1,'descend');
        if length(lmax)>N
            Imax=sort(Imax(1:N));
        else
            Imax=sort(Imax);
            N=length(lmax);
        end

        % We detect the lowest minima between two consecutive maxima
        bound=zeros(1,N);
        for i=1:N
           if i==1
               a=1;
           else
               a=Imax(i-1);
           end
           [lmin,ind]=sort(locmin(a:Imax(i)));
           tmp=lmin(1);
           n=1;
           if n<length(lmin)
               n=2;
               while ((n<=length(lmin)) && (tmp==lmin(n)))
                    n=n+1;
               end
           end
           bound(i)=a+ind(ceil(n/2))-2;
        end
    else
        k=1;
        for i=1:length(locmin)
           if locmin(i)<max(f2)
              bound(k) = i-1;
              k=k+1;
           end
        end
    end
    bound=[bound+band(1)-1,band(2)-1];
end


