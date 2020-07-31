%%
in2mm=25.4;
%Diameter = 1.5*in2mm+8; % mm
Diameter = 1.5*in2mm; % mm
Length = 2.9*in2mm; % mm
%sensor_id = 1:8;
L0=.414; L1=.989; L2=1.59; L3=2.162; L4=2.768;
Axial = Length-[L4-L0; L3-L0; L2-L0; L1-L0; L4-L0; L3-L0; L2-L0; L1-L0]*in2mm; % mm
Azimuth = [pi/4; pi/2; 3*pi/4; pi; 5*pi/4; 3*pi/2; 7*pi/4; 2*pi]; 
Axial_d = Length-[L4-L0; L3-L0; L2-L0; L1-L0; L4-L0; L3-L0; L2-L0; L1-L0]*in2mm; % mm
Azimuth_d = [3*pi/4; pi; 5*pi/4; 3*pi/2; 7*pi/4; 2*pi; pi/4; pi/2]; 
sensor_id = 8:-1:1;
channel_id = 0:7;
% sensor: 7, 8, 6, 5, 4, 3, 1, 2 
% channel: 0, 1, 3, 4, 5, 6, 7, 8
% plug_height = 10; % mm

%% plot cylinder
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

%plot sensor: channel 1-8
xr = Diameter/2 * sin(Azimuth);
yr = Diameter/2 * cos(Azimuth);
zr = Axial;
plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
for i = 1:8
    %text(xr(i)+ 2.5, yr(i), zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
    text(xr(i)+ 2.5, yr(i), zr(i),[num2str(channel_id(i)) '(' num2str(i+8) ')'],'fontsize',16)

end

% plot sensor: channel 1.-8.
xr_d = Diameter/2 * sin(Azimuth_d);
yr_d = Diameter/2 * cos(Azimuth_d);
zr_d = Axial_d;
plot3(xr_d, yr_d, zr_d, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
for i = 1:8
    %text(xr_d(i)+ 2.5, yr_d(i), zr_d(i),[num2str(channel_id(i)) '.(' num2str(sensor_id(i)) '.)'],'fontsize',16)
    text(xr_d(i)+ 2.5, yr_d(i), zr_d(i),[num2str(channel_id(i)) '.(' num2str(i) ')'],'fontsize',16)
end

% %% views
% subplot(2,2,1)
% R=[Diameter/2 Diameter/2];
% x0=0; y0=0; z0=0;
% N=100;
% [X,Y,Z] = cylinder(R,N);
% testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
% set(testsubjectO,'FaceAlpha',0.4,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
% axis equal
% set(gca,'fontname','arial','fontsize',16)
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% hold on
% 
% xr = Diameter/2 * sin(Azimuth);
% yr = Diameter/2 * cos(Azimuth);
% zr = Axial - 10;
% plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
% for i = 1:8
%     text(xr(i)+ 2.5, yr(i), zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
% end
% 
% subplot(2,2,2)
% R=[Diameter/2 Diameter/2];
% x0=0; y0=0; z0=0;
% N=100;
% [X,Y,Z] = cylinder(R,N);
% testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
% set(testsubjectO,'FaceAlpha',0.4,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
% axis equal
% set(gca,'fontname','arial','fontsize',16)
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% hold on
% xr = Diameter/2 * sin(Azimuth);
% yr = Diameter/2 * cos(Azimuth);
% zr = Axial - 10;
% plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
% for i = 1:8
%     text(xr(i)+ 2.5, yr(i), zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
% end
% view(0,0)
% 
% subplot(2,2,3)
% R=[Diameter/2 Diameter/2];
% x0=0; y0=0; z0=0;
% N=100;
% [X,Y,Z] = cylinder(R,N);
% testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
% set(testsubjectO,'FaceAlpha',0.4,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
% axis equal
% set(gca,'fontname','arial','fontsize',16)
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% hold on
% xr = Diameter/2 * sin(Azimuth);
% yr = Diameter/2 * cos(Azimuth);
% zr = Axial - 10;
% plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
% for i = 1:8
%     text(xr(i)+ 2.5, yr(i), zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
% end
% view(90,0)
% 
% subplot(2,2,4)
% R=[Diameter/2 Diameter/2];
% x0=0; y0=0; z0=0;
% N=100;
% [X,Y,Z] = cylinder(R,N);
% testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
% set(testsubjectO,'FaceAlpha',0.4,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
% axis equal
% set(gca,'fontname','arial','fontsize',16)
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% hold on
% r1=Diameter/2; %inch
% theta1=0:pi/18:2*pi;
% x1=r1*cos(theta1); y1=r1*sin(theta1);
% testsubject3 = patch(x1+x0,y1+y0,12*ones(size(x1)),'b');
% set(testsubject3,'FaceAlpha',0.4,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
% xr = Diameter/2 * sin(Azimuth);
% yr = Diameter/2 * cos(Azimuth);
% zr = Axial - 10;
% plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
% for i = 1:8
%     text(xr(i)+ 2.5, yr(i), zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
% end
% view(0,90)

%% views
figure
subplot(1,4,1)
R=[Diameter/2 Diameter/2];
x0=0; y0=0; z0=0;
N=100;
[X,Y,Z] = cylinder(R,N);
testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
set(testsubjectO,'FaceAlpha',0.2,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
axis equal
set(gca,'fontname','arial','fontsize',18)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
hold on

%plot sensor: channel 1-8
plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
for i = 1:8
    text(xr(i)+ 2.5, yr(i), zr(i),num2str(i+8),'fontsize',24)
end

% plot sensor: channel 1.-8.
plot3(xr_d, yr_d, zr_d, 'ko', 'markersize', 10, 'markerfacecolor','k','markeredgecolor','k')
for i = 1:8
    %text(xr_d(i)+ 2.5, yr_d(i), zr_d(i),[num2str(channel_id(i)) '.(' num2str(sensor_id(i)) '.)'],'fontsize',16)
    text(xr_d(i)+ 2.5, yr_d(i), zr_d(i),num2str(i),'fontsize',24)

end


subplot(1,4,2:3)
l = Diameter/2 * (2*pi - Azimuth);
l_d = Diameter/2 * (2*pi - Azimuth_d);
plot(l,zr,'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r'); hold on
plot(l_d,zr_d,'ko', 'markersize', 10, 'markerfacecolor','k','markeredgecolor','k')

axis equal
xlim([0 Diameter*pi])
ylim([0 Length])
testsubject3 = patch([0 Diameter*pi Diameter*pi 0],[0 0 Length Length],'b');
set(testsubject3,'FaceAlpha',0.2,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
set(gca,'fontname','arial','fontsize',18)
xlabel('L (mm)')
ylabel('Z (mm)')
for i = 1:8
    %text(l(i) + 2.5, zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
    text(l(i) + 2.5, zr(i),num2str(i+8),'fontsize',24)

end
for i = 1:8
    %text(l_d(i) + 2.5, zr_d(i),[num2str(channel_id(i)) '.(' num2str(sensor_id(i)) '.)'],'fontsize',16)
    text(l_d(i) + 2.5, zr_d(i),num2str(i),'fontsize',24)
end
grid on

subplot(1,4,4)
R=[Diameter/2 Diameter/2];
x0=0; y0=0; z0=0;
N=100;
[X,Y,Z] = cylinder(R,N);
testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
set(testsubjectO,'FaceAlpha',0.4,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
axis equal
set(gca,'fontname','arial','fontsize',18)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
hold on
r1=Diameter/2; %inch
theta1=0:pi/18:2*pi;
x1=r1*cos(theta1); y1=r1*sin(theta1);
testsubject4 = patch(x1+x0,y1+y0,Length*ones(size(x1)),'b');
set(testsubject4,'FaceAlpha',0.2,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
% plot sensor: channel 0-7
plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
for i = 1:8
    %text(xr(i)+ 2.5, yr(i)+2.5, zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
    text(xr(i)+ 2.5, yr(i)+5, zr(i),num2str(i+8),'fontsize',24)
end

% plot sensor: channel 0.-7.
plot3(xr_d, yr_d, zr_d, 'ko', 'markersize', 10, 'markerfacecolor','k','markeredgecolor','k')
for i = 1:8
    %text(xr_d(i)+ 2.5, yr_d(i), zr_d(i),[num2str(channel_id(i)) '.(' num2str(sensor_id(i)) '.)'],'fontsize',16)
    text(xr_d(i)+ 5, yr_d(i), zr_d(i),num2str(i),'fontsize',24)
end

view(0,90)
box on