%% 
% G1: (-6.7159, -6.5545) mm
% G2: ( 0.7749, -9.3313) mm
% G3: (-7.7814,  5.0047) mm
% G4: ( 5.0369,  8.1366) mm
% G5: ( 8.8792, -1.1947) mm

%% plot cylinder
in2mm=25.4;
Diameter = 1.5*in2mm+8; % mm
Length = 3*in2mm; % mm
%sensor_id = 1:8;
Axial = [12/5; 9/5; 6/5; 3/5; 12/5; 9/5; 6/5; 3/5]*in2mm; % mm
Azimuth = [pi/2; pi/4; 0; 7*pi/4; 3*pi/2; 5*pi/4; pi; 3*pi/4];
Axial_d = [12/5; 9/5; 6/5; 3/5; 12/5; 9/5; 6/5; 3/5]*in2mm; % mm
Azimuth_d = [0; 7*pi/4; 3*pi/2; 5*pi/4; pi; 3*pi/4; pi/2; pi/4];
sensor_id = 8:-1:1;
channel_id = 0:7;

Gx=[-6.7159; 0.7749; -7.7814; 5.0369; 8.8792];
Gy=[-6.5545; -9.3313; 5.0047; 8.1366; -1.1947];
Gz=Length*ones(5,1);

figure
R=[Diameter/2 Diameter/2];
x0=0; y0=0; z0=0;
N=100;
[X,Y,Z] = cylinder(R,N);
testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
set(testsubjectO,'FaceAlpha',0.2,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
axis equal
set(gca,'fontname','arial','fontsize',16)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
hold on
hold on
r1=Diameter/2; %inch
theta1=0:pi/18:2*pi;
x1=r1*cos(theta1); y1=r1*sin(theta1);
testsubject3 = patch(x1+x0,y1+y0,0*ones(size(x1)),'b');
set(testsubject3,'FaceAlpha',0.2,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)

color={'r','g','b','k','m'};
for Ge=1:5
    plot3(Gx(Ge),Gy(Ge),Gz(Ge), 'd', 'markersize', 14, 'markerfacecolor',color{Ge},'markeredgecolor',color{Ge})
end

%% First ev of each group
flow=10e3;
fhigh=300e3;
fs=25e6;
[a,b] = butter(2,[flow*2/fs fhigh*2/fs],'bandpass');
t=0:1/fs:1/fs*2500;
for G=10
    figure;
    sig=filtfilt(a,b,ev_wfm{G});
    for ch=0:7  
        subplot(16,1,ch*2+1:ch*2+2)
        plot(t*1e6,sig(:,ch+1)); hold on
        plot([tpt{G}(ch+1,1) tpt{G}(ch+1,1)]*1e6, [-max(abs(sig(:,ch+1))) max(abs(sig(:,ch+1)))],'k')
        if mod(G,2)==1
            theoP=t0_all{G}(1)+sqrt((Gx(ceil(G/2))-xr(ch+1))^2+(Gy(ceil(G/2))-yr(ch+1))^2+(Gz(ceil(G/2))-zr(ch+1))^2)/Vp/1e3;
        else
            theoP=t0_all{G}(1)+sqrt((Gx(ceil(G/2))-xr_d(ch+1))^2+(Gy(ceil(G/2))-yr_d(ch+1))^2+(Gz(ceil(G/2))-zr_d(ch+1))^2)/Vp/1e3;
        end
        plot([theoP theoP]*1e6, [-max(abs(sig(:,ch+1))) max(abs(sig(:,ch+1)))],'r')
        xlim([t(1) t(end)]*1e6)
        if ch==7
            set(gca,'fontsize',14,'fontname','arial')
            xlabel('Time (\musec)')
        else
            set(gca,'fontsize',14,'fontname','arial','xticklabel',[])
        end
        if mod(G,2)==0
            ylabel(['ch ' num2str(ch) '.'])
        else
            ylabel(['ch ' num2str(ch)])
        end
    end
end