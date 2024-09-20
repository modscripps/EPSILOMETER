function f=epsilon2_correct(eps_res,kvis,kmin,kmax)

pf=[-3.12067617e-05, -1.31492870e-03, -2.29864661e-02, -2.18323426e-01, -1.23597906, ...
    -4.29137352,-8.91987933, -9.58856889, -2.41486526];

ler=log10(eps_res);
	
if ler <= -5 % adjust shear with missing variance
	lf_eps=0;

    % % [kpan,Ppan] = panchev(eps_res,kvis);
    % [Ppan,kpan]=nasmyth(eps_res,kvis);
    % idx_kmax2=find(kpan>=kmax,1,'first');
    % idx_kmax1=find(kpan>=kmin,1,'first');
    % idx_range=idx_kmax1:idx_kmax2;
    % tot_var=trapz(kpan(idx_kmax1:end),Ppan(idx_kmax1:end));
    % 
    % Pan_var_explain=trapz(kpan(idx_range),Ppan(idx_range));
    % lf_eps=log10(2-Pan_var_explain/tot_var); 

elseif ler>-5 && ler<=-1 % range of fit
	lf_eps=polyval(pf,ler);	
	if ler < -2 % apply viscosity correction
 		lf_eps=lf_eps+0.05*(ler+6)*(1.3e-6-kvis)/(0.3e-6);
	end
elseif ler>-1
	lf_eps=polyval(pf,-1);
end

f=10.^lf_eps;

