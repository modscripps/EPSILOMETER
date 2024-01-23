function vnav=vnav_pitch_roll_heading(vnav)

    MagZ=vnav.compass(:,3);
    MagY=vnav.compass(:,2);
    MagX=vnav.compass(:,1);

    % Now if heading ~ 0; X direction aligned with North
    Avx=vnav.acceleration(:,1); % m/s^-2
    Avy=vnav.acceleration(:,2); % m/s^-2
    Avz=vnav.acceleration(:,3); % m/s^-2
    
    %convert to XYZ
    pitch = atan2 ( Avx, sqrt(Avy.^2 + Avz.^2));
%     roll  = atan2 ( -Avy , sqrt(Avx.^2 + Avz.^2));
    roll  = atan2 ( -Avy , Avz);

    X=MagZ.*sin(roll) - MagY.*cos(roll);
    Y=MagX.*cos(pitch)+MagY.*sin(pitch).*sin(roll)+MagZ.*sin(pitch).*cos(roll);
    yaw=atan2(X,Y);
    heading=yaw;

    Av=[Avx Avy Avz].';
    %theta=pitch psi=yaw phi=roll
    DCM(1,1,:)=  cos(pitch).*cos(yaw);
    DCM(1,2,:)=  cos(pitch).*sin(yaw);
    DCM(1,3,:)= -sin(pitch);

    DCM(2,1,:)=  sin(roll).*sin(pitch).*cos(yaw) - cos(roll).*sin(yaw);
    DCM(2,2,:)=  sin(roll).*sin(pitch).*sin(yaw) + cos(roll).*cos(yaw);
    DCM(2,3,:)=  sin(roll).*cos(pitch);

    DCM(3,1,:)=  cos(roll).*sin(pitch).*cos(yaw) + sin(roll).*sin(yaw);
    DCM(3,2,:)=  cos(roll).*sin(pitch).*sin(yaw) - sin(roll).*cos(yaw);
    DCM(3,3,:)=  cos(roll).*cos(pitch);

    ANED=zeros(3,length(pitch));
    for t=1:length(pitch)
        R=squeeze(DCM(:,:,t));
        ANED(:,t)=R*Av(:,t);
    end

    vnav.pitch    = pitch;
    vnav.roll     = roll;
    vnav.heading  = heading;
    vnav.ANED      = ANED; % Acceleration North East Down
    
    
    

