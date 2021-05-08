function ax1=add_scan_to_plot_epsi_chi_spectra(scan,scan_name,ax1)

for i=1:length(ax1)
    hold(ax1(i),'on');
end


%%

a=3;
loglog(ax1(a),scan.fe,scan.Pt,'m','linewidth',2)
ax1(a).Legend.String={ax1(a).Legend.String{1:end-1}, ['Pt ' scan_name]};
ax1(a).XLim=[.3 160];
ax1(a).YLim=[1e-10 1e-4];
ax1(a).FontSize=20;


a=1;
loglog(ax1(a),scan.fe,scan.Pv,'m','linewidth',2)
ax1(a).Legend.String={ax1(a).Legend.String{1:end-1}, ['Pv ' scan_name]};
ax1(a).XLim=[.3 160];
ax1(a).YLim=[1e-11 1e-6];
ax1(a).FontSize=20;


a=4;
loglog(ax1(a),scan.fe,scan.Pa2,'g','linewidth',2)
loglog(ax1(a),scan.fe,scan.Pa3,'m','linewidth',2)
ax1(a).Legend.String={ax1(a).Legend.String{1:end-2}, ['Pa2 ' scan_name],['Pa3 ' scan_name]};
ax1(a).XLim=[.3 160];
ax1(a).YLim=[1e-8 1e-4];
ax1(a).FontSize=20;
ax1(a).FontName='Times New Roman';


a=2;
semilogx(ax1(a),scan.fe,scan.Cu1a3,'m','linewidth',2)
ax1(a).Legend.String={ax1(a).Legend.String{1:end-1}, ['s-a3 ' scan_name]};
grid(ax1(a),'on')
ax1(a).XLim=[.3 160];
ax1(a).FontSize=20;


for i=1:length(ax1)
    set(ax1(i),'children',flipud(get(ax1(i),'children')))
    hold(ax1(i),'off');
end


