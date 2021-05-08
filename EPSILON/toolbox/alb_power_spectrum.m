function P11=alb_power_spectrum(data,dk)
    try 
        [nbscan,Lscan]=size(data);
        if nbscan>Lscan
            warning('data dimension should be data(segment,time)')
        end
        window = ones(nbscan,1)*hanning(Lscan).';
        wc2=1/mean(window(1,:).^2);            % window correction factor
        %datap  = window.*(data- mean(data,2)* ones(1,Lscan));
        datap=window.* detrend(data.').';
        P1  = fft(datap,[],2);
        P11 = conj(P1).*P1./Lscan^2/dk*wc2;
    catch 
        disp('data should be 2 dimension (nb segment, time)')
    end
