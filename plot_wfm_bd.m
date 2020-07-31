%%
fs=12.5e6; % 12.5 MHz
f0=100e3; % 100 k
flow = 50e3;     % min freq of passband, Hz
fhigh = 300e3;    % max freq of passband, Hz
[b,a]=butter(2,[flow*2/fs fhigh*2/fs],'bandpass');
wfm=zeros(4001,16);
t=0:1/fs:4000/fs;
for G=1:2
    for st=1:8
        u=filtfilt(b,a,ev_wfm{G}(:,st));
        u=u(200:4200);  
        wfm(:,8*(G-1)+st)=u;
    end  
end
%%
max_wfm=max(max(abs(wfm(:,9:12))));
figure;
for i=1:4
    plot(t*1e6,wfm(:,8+i)/max_wfm+i*2,'b','linewidth',2); hold on
end
set(gca,'fontsize',18,'ytick',2:2:8,'yticklabel',{'ch 0','ch 1','ch 2','ch 3'})
xlabel('Time (\mus)')
xlim([0 300])