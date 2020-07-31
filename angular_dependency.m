%% sensor radiation
theta = 0:90;
a=5e-3; %m
vp=6011.6; % m/s
vs=3093.0; % m/s
sigma_zz=1; R = 1;
h1=figure;
for f=100e3:300e3:1000e3 % hz
    c=[f/1000e3 f/1000e3 f/1000e3]*.6;
    k_alpha = 2*pi*f/vp;
    j1_over_x=besselj(1,k_alpha*a*sin(deg2rad(theta)))./(k_alpha*a*sin(deg2rad(theta)));
    j1_over_x(1)=1/2;
    u_R = a^2*sigma_zz/4/pi*exp(-1i*k_alpha*R)/R*2*j1_over_x.*...
        ((vs/vp)^2*cos(deg2rad(theta)).*(1-2*(vs/vp)^2*sin(deg2rad(theta)).^2))./...
        ((1-2*(vs/vp)^2*sin(deg2rad(theta)).^2).^2 + 4*(vs/vp)^3*sin(deg2rad(theta)).^2.*cos(deg2rad(theta)).*sqrt(1-(vs/vp)^2*sin(deg2rad(theta)).^2));
    u_R_normlize=abs(u_R)/max(abs(u_R));
    plot(u_R_normlize.*sin(deg2rad(theta)),u_R_normlize.*cos(deg2rad(theta)),'color',c,'linewidth',2);
    hold on
    plot(u_R_normlize.*sin(deg2rad(-theta)),u_R_normlize.*cos(deg2rad(-theta)),'color',c,'linewidth',2);
end
plot(cos(deg2rad(theta)).*sin(deg2rad(theta)),cos(deg2rad(theta)).*cos(deg2rad(theta)),'r','linewidth',2)
hold on
plot(cos(deg2rad(theta)).*sin(deg2rad(-theta)),cos(deg2rad(theta)).*cos(deg2rad(-theta)),'r','linewidth',2);

axis off
%%
h2=figure;
for f=100e3:300e3:1000e3 % hz
    c=[f/1000e3 f/1000e3 f/1000e3]*.6;
    k_alpha = 2*pi*f/vp;
    j1_over_x=besselj(1,k_alpha*a*sin(deg2rad(theta)))./(k_alpha*a*sin(deg2rad(theta)));
    j1_over_x(1)=1/2;
    u_R = a^2*sigma_zz/4/pi*exp(-1i*k_alpha*R)/R*2*j1_over_x.*...
        ((vs/vp)^2*cos(deg2rad(theta)).*(1-2*(vs/vp)^2*sin(deg2rad(theta)).^2))./...
        ((1-2*(vs/vp)^2*sin(deg2rad(theta)).^2).^2 + 4*(vs/vp)^3*sin(deg2rad(theta)).^2.*cos(deg2rad(theta)).*sqrt(1-(vs/vp)^2*sin(deg2rad(theta)).^2));
    u_R_normlize=abs(u_R)/max(abs(u_R));
    plot(theta,u_R_normlize,'color',c,'linewidth',2); hold on
end
plot(theta,cos(deg2rad(theta)),'r','linewidth',2)
set(gca,'fontsize',18,'xtick',0:15:90)
xlabel('Incident angle (degree)')
ylabel('Normalized amplitude')
legend('100 kHz','400 kHz', '700 kHz', '1 MHz','cosine function')