function PlotYEIonFCTD2(FCTD)
disp(date);
disp(FCTD);
disp(datestr(FCTD.time(end)));
% date: 2013 05 07
pts = 96*2;

if pts > numel(FCTD.time)
    pts = numel(FCTD.time);
end
applyFiltering = false;
if applyFiltering
    comp = FCTD.compass(end-fpts-1:end,:);
    gyro = FCTD.gyro(end-fpts-1:end,:);
    acceleration = FCTD.acceleration(end-fpts-1:end,:);
    
    FCTD.compass = medfilt1(FCTD.compass,fpts/2,[],1);
    FCTD.gyro = medfilt1(FCTD.gyro,fpts/2,[],1);
    FCTD.acceleration = medfilt1(FCTD.acceleration,fpts/2,[],1);
    
    mygausswin = gausswin(fpts);
    mygausswin = mygausswin/sum(mygausswin);
    
    FCTD.compass = conv2(FCTD.compass,mygausswin,'same');
    FCTD.gyro = conv2(FCTD.gyro,mygausswin,'same');
    FCTD.acceleration = conv2(FCTD.acceleration,mygausswin,'same');
    
    FCTD.compass(end-fpts-1:end,:) = comp;
    FCTD.gyro(end-fpts-1:end,:) = gyro;
    FCTD.acceleration(end-fpts-1:end,:) = acceleration;
    
end

figure(1000);
set(gcf,'renderer','painters');
% clf;
% % subplot(2,3,5:6);
% % plot(FCTD.time(i+(0:(pts-1))),FCTD.pressure(i+(0:(pts-1))),'r-*','markersize',2);
% % ylim([0 2000]);
% axis ij;
% % hold on;
% ylabel('Pressure [dbar]','interpreter','latex');
% xlabel('Time (MM:SS) [UTC]','interpreter','latex');
% grid on;
% box on;
% datetick('x','MM:SS','keeplimits')
LHS = true;
if LHS
    multiplier = -1;
else
    multiplier = 1;
end

[Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis([0, 1,0]);
Rot_Mat = @(p,t,s)[ cos(t)*cos(s), -cos(p)*sin(s) + sin(p)*sin(t)*cos(s),  sin(p)*sin(s) + cos(p)*sin(t)*cos(s);
    cos(t)*sin(s),  cos(p)*cos(s) + sin(p)*sin(t)*sin(s), -sin(p)*cos(s) + cos(p)*sin(t)*sin(s);
    -sin(t),         sin(p)*cos(t),                         cos(p)*cos(t)];
acceleration = (Rot_Mat(0,0,pi/2)*(Rot_mat*(FCTD.acceleration(end-pts+1:end,:)')))';
acceleration(:,3) = -multiplier*acceleration(:,3);

acc_xy_length = sqrt(sum(acceleration(:,1:2).^2,2));
acc_length = sqrt(sum(acceleration.^2,2));

%pitch
subplot(2,3,1);
plot([0, 0],[-1, 1],'color','k','linewidth',1,'linestyle','--');
hold on;
plot([-1, 1],[0, 0],'color','k','linewidth',1,'linestyle','--');
plot([-1, 1],[1, -1],'color','k','linewidth',1,'linestyle','--');
plot([1, -1],[1, -1],'color','k','linewidth',1,'linestyle','--');

arrow_length = sqrt(acceleration(:,1).^2+acceleration(:,3).^2);
plot(acceleration(:,3)./arrow_length,...
    acceleration(:,1)./arrow_length,...
    'color','b','linewidth',0.5);
quiver(0,0,acceleration(end,3)./arrow_length(end),...
    acceleration(end,1)./arrow_length(end),0,...
    'color','b','linewidth',5);

text(0,0,'\qquad {\bf horizontal} \qquad',...
    'backgroundcolor','none','color','k','interpreter','latex','fontweight','bold',...
    'horizontalalignment','center','verticalalignment','middle','rotation',0,...
    'fontsize',15);
text(0,-1.45,sprintf('\\quad%0.0f$^\\circ$ with respect to horizontal\\qquad',...
    atan2(acceleration(end,1),acceleration(end,3))*180/pi),...
    'backgroundcolor','k','color','w','interpreter','latex','fontweight','bold',...
    'horizontalalignment','center','verticalalignment','middle');

xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('$x$','interpreter','latex');
ylabel('$z$','interpreter','latex');
set(gca,'dataaspectratio',[1 1 1],'xtick',[],'ytick',[]);
grid on;
box on;
hold off;
title('Pitch of FCTD Fish [$r$-$z$ plane]','interpreter','latex');

% roll
subplot(2,3,2);
plot([0, 0],[-1, 1],'color','k','linewidth',1,'linestyle','--');
hold on;
plot([-1, 1],[0, 0],'color','k','linewidth',1,'linestyle','--');
plot([-1, 1],[1, -1],'color','k','linewidth',1,'linestyle','--');
plot([1, -1],[1, -1],'color','k','linewidth',1,'linestyle','--');

text(0,0,'\qquad {\bf vertical} \qquad',...
    'backgroundcolor','none','color','k',...
    'interpreter','latex','fontweight','bold',...
    'horizontalalignment','center',...
    'verticalalignment','middle','rotation',90,...
    'fontsize',15);
arrow_length = sqrt(acceleration(:,2).^2+acceleration(:,3).^2);
plot(acceleration(:,2)./arrow_length,...
    acceleration(:,3)./arrow_length,...
    'color','b','linewidth',0.5);
quiver(0,0,acceleration(end,2)./arrow_length(end),...
    acceleration(end,3)./arrow_length(end),0,...
    'color','b','linewidth',5);
hold on;
text(0,-1.45,sprintf('\\quad%0.0f$^\\circ$ with respect to vertical\\qquad',...
    atan2(acceleration(end,2),acceleration(end,3))*180/pi),...
    'backgroundcolor','k','color','w','interpreter','latex','fontweight','bold',...
    'horizontalalignment','center','verticalalignment','middle');

xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('$x$','interpreter','latex');
ylabel('$z$','interpreter','latex');
set(gca,'dataaspectratio',[1 1 1],'xtick',[],'ytick',[]);
grid on;
box on;
hold off;
title('Roll of FCTD Fish [$y$-$z$ plane]','interpreter','latex');


comp = FCTD.compass(end-pts+1:end,:);
gyro = FCTD.gyro(end-pts+1:end,:);
acceleration = FCTD.acceleration(end-pts+1:end,:);

comp(:,3) = multiplier*comp(:,3);
gyro(:,3) = multiplier*gyro(:,3);
acceleration(:,3) = multiplier*acceleration(:,3);


for k = 1:pts
    [Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis(acceleration(k,:));
    comp(k,:) = Rot_mat*(comp(k,:)');
    gyro(k,:) = Rot_mat*(gyro(k,:)');
    acceleration(k,:) = Rot_mat*(acceleration(k,:)');
end
% [Phi, Theta, Psi, Rot_mat] = SN_RotateToZAxis([0, 1,0]);
% Rot_Mat = @(p,t,s)[ cos(t)*cos(s), -cos(p)*sin(s) + sin(p)*sin(t)*cos(s),  sin(p)*sin(s) + cos(p)*sin(t)*cos(s);
%     cos(t)*sin(s),  cos(p)*cos(s) + sin(p)*sin(t)*sin(s), -sin(p)*cos(s) + cos(p)*sin(t)*sin(s);
%     -sin(t),         sin(p)*cos(t),                         cos(p)*cos(t)];
% comp = (Rot_Mat(0,0,pi/2)*(Rot_mat*(comp')))';
% gyro = (Rot_Mat(0,0,pi/2)*(Rot_mat*(gyro')))';
% acceleration = (Rot_Mat(0,0,pi/2)*(Rot_mat*(acceleration')))';
if(pts>24)
    comp_2_plot1 = nanmedfilt1(comp,12,[],1);
end
Comp_2_plot1 = repmat(sqrt(sum(comp_2_plot1(:,1:2).*comp_2_plot1(:,1:2),2)),[1 3]);
comp_2_plot1 = comp_2_plot1./Comp_2_plot1;



subplot(2,3,4);
cla;
h = polar([0 2*pi], [0 1.25]);
delete(h)
hold on;
% keyboard;
% h = polar(atan2(comp(:,2),comp(:,1)),(1:size(comp,1))'/size(comp,1));
h = polar(atan2(comp_2_plot1(:,2),comp_2_plot1(:,1)),(1:size(comp_2_plot1,1))'/size(comp_2_plot1,1));
% keyboard;
hold on;
h1 = compass(comp_2_plot1(end,1),comp_2_plot1(end,2),'b');
hold off
set(h,'color','b','marker','.','markersize',2,'markeredgecolor','b','linestyle',':');
set(h1,'color','r','linewidth',5);
set(gca,'yTicklabel',{});
get(gca,'yTicklabel')


% acceleration time plot
subplot(2 ,3,[5 6]) 
plot(FCTD.time(end-pts+1:end),acceleration./repmat(sqrt(sum(acceleration.^2,2)),[1 3]),'linewidth',2)
hold on;
plot(FCTD.time(end-pts+1:end),...
    comp./repmat(sqrt(sum(comp.^2,2)),[1 3]),...
    '--o','linewidth',2)
hold off;
datetick('x','keeplimits');
legend('acc_x','acc_y','acc_z','mag_x','mag_y','mag_z');

subplot(2,3,3);
cla;
grid off;
box off;
set(gca,'xtick',[],'ytick',[],'xcolor','w','ycolor','w');
ylim([0 1]);
xlim([0,1]);
axis ij;
text(0,0,'Top panels: pitch and roll of FCTD Fish','interpreter','latex');
hold on;
text(0,0.075,'Bottom panels: acc. is projected to $z$-axis','interpreter','latex');
text(0,0.15,'[blue - COMPASS, red - GYRO, green - ACC]','interpreter','latex');
text(0,0.255,'Use roll of fish to figure out twists on cable','interpreter','latex');
text(0,0.33,'when pitch is less than 75$^\circ$ from vertical,','interpreter','latex');
text(0,0.405,'otherwise, use bottom-left panel.','interpreter','latex');

text(0,0.705,'Current CTD data:','interpreter','latex');
text(0,0.810,...
    sprintf('Pressure:%4.01f, Temperature: %2.02f, Conductivity: %2.02f',...
    FCTD.pressure(end),FCTD.temperature(end),FCTD.conductivity(end)),...
    'interpreter','latex');

text(0,0.90,['Oldest record: ' datestr(FCTD.time(1),'yyyy-mm-dd HH:MM:SS.FFF')],'interpreter','latex','fontsize',12);
text(0,0.95,['Newest record: ' datestr(FCTD.time(end),'yyyy-mm-dd HH:MM:SS.FFF')],'interpreter','latex','fontsize',12);
text(0,1,['Last update: ' datestr(now,'yyyy-mm-dd HH:MM:SS.FFF')],'interpreter','latex','fontsize',12);
hold off;
grid on;
box on;
end