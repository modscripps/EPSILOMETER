function M = read_digilient_network_analysis(filename)
    M1 = csvread(filename,6);
    M.f=M1(:,1); % Hz
    M.Phase = M1(:,4);  %Deg
    M.Vin   = 10.^(M1(:,2)/20);  %Volt
    M.Vout  = 10.^(M1(:,3)/20);  %Volt
    
    % dB=20 log_10 (Vout/Vin) 