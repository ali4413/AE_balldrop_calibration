%%
[spec_bd,f]=plotfft(wfm,1/fs);
clr={'b','r','k','g'};
for i=1:4
    semilogx(f/1e3,20*log10(abs(spec_bd(:,8+i))/max(max(abs(spec_bd(:,9:12))))),'color',clr{i},'linewidth',2); hold on
    xlim([f(1) f(end)]/1e3)
    set(gca,'fontsize',16)
    xlabel('Frequency (kHz)')
    ylabel('Spectra Amplitude (dB)')
end
legend('ch 0', 'ch 1', 'ch 2', 'ch 3')
grid on

%% spectra of force pulse
mu1=0.3; E1=200e9; delta1=(1-mu1^2)/(pi*E1);
mu2=0.35; E2=100e9; delta2=(1-mu2^2)/(pi*E2);
rho1=4.51e3; % kg/m^3
m1=1.0625e-3; % kg
R1=(m1/rho1*3/4/pi)^(1/3); % 6.36 mm; tunnel diameter 7.00 mm; hight 110.40 mm
v0=sqrt(2*9.8*.11); %m/s
tc=4.53*(4*rho1*pi*(delta1+delta2)/3)^(2/5)*R1*v0^(-1/5);
Fmax=1.917*rho1^(3/5)*(delta1+delta2)^(-2/5)*R1^2*v0^(6/5);
t_load=linspace(0,tc,100);
f_load=Fmax*abs((sin(pi*t_load'/tc)).^(3/2));
figure; plot(t_load*1e6,f_load,'b','linewidth',2)
set(gca,'fontsize',18)
xlabel('Time (\mus)')
ylabel('Force (N)')
%%
[spec_source,ffc]=plotfft([f_load; zeros(2049-length(f_load),1)],1/fs);
semilogx(f/1e3,20*log10(abs(spec_source)/max(max(abs(spec_source)))),'color','r','linewidth',2); hold on
xlim([f(1) f(end)]/1e3)
set(gca,'fontsize',16)
xlabel('Frequency (kHz)')
ylabel('Spectra Amplitude (dB)')

%%
TF_bd=zeros(2049,16);
waterlevel = 0;
for i=1:16
    TF_bd(:,i)=sqrt((spec_bd(:,i).*conj(spec_source))./((abs(spec_source)).^2 + waterlevel.*repmat(mean(abs(spec_source).^2),size(spec_source,1),1)));
end
%%
figure
for i=1:4
    semilogx(f/1e3,20*log10(abs(TF_bd(:,8+i))/max(max(abs(TF_bd(:,9:12))))),'-o','color',clr{i},'linewidth',2); hold on
    xlim([f(1) f(end)]/1e3)
    set(gca,'fontsize',16)
    xlabel('Frequency (kHz)')
    ylabel('Spectra Amplitude (dB)')
end
legend('ch 0', 'ch 1', 'ch 2', 'ch 3')
grid on
title('Instrument-apparatus response spectra')