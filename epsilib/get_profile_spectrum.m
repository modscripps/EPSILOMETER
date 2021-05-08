function [k,P1,P11,Co12]=get_profile_spectrum(data,k)
%
%  input: data
% . data : epsi data
% . k:  frequency array 
%  Created by Arnaud Le Boyer on 7/28/18.


    switch length(size(data))
        case 3 % reshape into 2D matrice for fft and then reshape
            [nb_sensor,nb_scan,Lscan]=size(data);
            data=reshape(data,[nb_sensor* nb_scan Lscan]);
            Lax1=nb_sensor* nb_scan;
            size_data=3;
        case 2
            [Lax1,Lscan]=size(data);
            size_data=2;
        otherwise
            warning('no valid size for data : get power spectrum')
    end
    
    dk=k(1);
    window = ones(Lax1,1)*hanning(Lscan).';
    wc2=1/mean(window(1,:).^2);            % window correction factor
    datap  = window.*(data- mean(data,2)* ones(1,Lscan));
    P1  = fft(datap,[],2);
    P11 = conj(P1).*P1./Lscan^2/dk*wc2;
    
    if size_data==3
        P1=reshape(P1,[nb_sensor,nb_scan,Lscan]);
        P11=reshape(P11,[nb_sensor,nb_scan,Lscan]);
        P12=zeros(nb_sensor,nb_sensor-1,nb_scan,Lscan);
        Co12=zeros(nb_sensor,nb_sensor-1,nb_scan,Lscan);
        ind_nbsensor=1:nb_sensor;
        for j=1:nb_sensor
            tempo=shiftdim(repmat(squeeze(P1(j,:,:)),[1,1,nb_sensor-1]),2);
            P12(j,:,:,:)=conj(tempo).*P1(ind_nbsensor~=j,:,:)./Lscan^2/dk*wc2;
            tempo=shiftdim(repmat(squeeze(P11(j,:,:)),[1,1,nb_sensor-1]),2);
            Co12(j,:,:,:)=squeeze(P12(j,:,:,:)).^2./(tempo.*P11(ind_nbsensor~=j,:,:));
        end
    end
    
    if rem(Lscan,2)==0
        k=-Lscan/2*dk:dk:Lscan/2*dk-dk;
    else
        kp=dk:dk:dk*floor(Lscan/2);
        k=[fliplr(-kp) 0 kp];
    end
    k=fftshift(k);
