function plot_binned_turbulence(Meta_Data,MS)

% if MS is not an input
if nargin==1
    load(fullfile(Meta_Data.L1path,'Turbulence_Profiles.mat'),'MS')
end

% compite binned epsilon for all profiles
Epsilon_class=calc_binned_epsi(MS);
Chi_class=calc_binned_chi(MS);

% plot binned epsilon for all profiles
close all
F1=figure(1);F2=figure(2);F3=figure(3);F4=figure(4);F5=figure(5);F6=figure(6);
[F1,F2]=plot_binned_epsilon(Epsilon_class,[Meta_Data.mission '-' Meta_Data.deployment],F1,F2,Meta_Data);
[F3,F4]=plot_binned_epsilon_co(Epsilon_class,[Meta_Data.mission '-' Meta_Data.deployment ' coh'],F3,F4,Meta_Data);
% [F5,F6]=plot_binned_epsilon_eof(Epsilon_class,[Meta_Data.mission '-' Meta_Data.deployment ' eof'],F5,F6,Meta_Data);
F1.PaperPosition = [0 0 30 20];F2.PaperPosition = [0 0 30 20];
F3.PaperPosition = [0 0 30 20];F4.PaperPosition = [0 0 30 20];
% F5.PaperPosition = [0 0 30 20];F6.PaperPosition = [0 0 30 20];
print(F1,fullfile(Meta_Data.L1path,[Meta_Data.deployment '_binned_epsilon1_t3s.png']),'-dpng2')
print(F2,fullfile(Meta_Data.L1path,[Meta_Data.deployment '_binned_epsilon2_t3s.png']),'-dpng2')
print(F3,fullfile(Meta_Data.L1path,[Meta_Data.deployment '_binned_epsilon1_co_t3s.png']),'-dpng2')
print(F4,fullfile(Meta_Data.L1path,[Meta_Data.deployment '_binned_epsilon2_co_t3s.png']),'-dpng2')
% print(F5,fullfile(Meta_Data.L1path,[Meta_Data.deployment '_binned_epsilon1_eof_t3s.png']),'-dpng2')
% print(F6,fullfile(Meta_Data.L1path,[Meta_Data.deployment '_binned_epsilon2_eof_t3s.png']),'-dpng2')


%% Chi plot
switch Meta_Data.MAP.temperature
    case 'Tdiff'
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_noise.mat'),'n0','n1','n2','n3');
    otherwise
        FPO7noise=load(fullfile(Meta_Data.CALIpath,'FPO7_notdiffnoise.mat'),'n0','n1','n2','n3');
end
dTdV(1)=Meta_Data.epsi.t1.dTdV; %1/0.025 V/deg 
dTdV(2)=Meta_Data.epsi.t2.dTdV; %1/0.025 V/deg 

f=1/3:1/3:320/2;
logf=log10(f);
w_th=.6;
h_freq=get_filters_MADRE(Meta_Data,f);
noise=(2*pi*f./w_th).^2.*10.^(FPO7noise.n0+FPO7noise.n1.*logf+ ...
           FPO7noise.n2.*logf.^2+...
           FPO7noise.n3.*logf.^3).*dTdV(1).^2./h_freq.FPO7(w_th).*w_th;

%i=1:10;j=2;
figure(7)
i=1:12;j=12;
l1=loglog(Chi_class.k,squeeze(Chi_class.Pbatch21(i,j,:)),'Color',.8*[1 1 1],'linewidth',2);
hold on;
l2=loglog(Chi_class.k,squeeze(Chi_class.mPphiT21(i,j,:)),'linewidth',2);
l3=loglog(f./w_th,3.*noise,'Color',[.2 .2 .2],'linewidth',2);
grid on
xlabel('k [cpm]')
ylabel('$(\phi^T_k)^2$  [ $(^{\circ}C . m^{-1})^2$ / cpm]','interpreter','latex')
set(gca,'fontsize',20)
legend([l1(1) l2(1) l3],{'batchelor','data','noise'},'location','best')
print(fullfile(Meta_Data.L1path,[Meta_Data.deployment '_j12_binned_chi_increpsi_t3s.png']),'-dpng2')

figure(8)
i=9;j=1:12;
loglog(Chi_class.k,squeeze(Chi_class.Pbatch21(i,j,:)),'Color',.8*[1 1 1],'linewidth',2)
hold on;
loglog(Chi_class.k,squeeze(Chi_class.mPphiT21(i,j,:)),'linewidth',2);
grid on
xlabel('k [cpm]')
ylabel('$(\phi^T_k)^2$  [ $(^{\circ}C . m^{-1})^2$ / cpm]','interpreter','latex')
set(gca,'fontsize',20)
print(fullfile(Meta_Data.L1path,[Meta_Data.deployment '_binned_chi_incrchi_t3s.png']),'-dpng2')

