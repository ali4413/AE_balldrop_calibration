DirM = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc/mat_bd/bd_cp30_ds10';
DirS = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc/wfm_plot_continous';
DirD = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc';
fs = 12.5e6;
bd_all = [repmat(6,1,8) repmat(8,1,8)];
ch_all = repmat(0:7,1,2);
dt=1/fs;
fn=1/(2*dt);
f1=5e3; f2=500e3;
[af,bf]=butter(2,[f1/fn  f2/fn]);
%t_detect = cell(16,1);
offset=3.15;
for test = 3
    load([DirD '/bd_cp30_ds10_' num2str(test) '.tdms_stalta.mat'])
    for ch_idx = 1:16
        h=figure;
        for s_idx = 5
            filename = ['bd_cp30_ds10_' num2str(test) '.tdms_PXI1Slot' num2str(bd_all(ch_idx)) '_' num2str(s_idx) '.mat'];
            load([DirM '/' filename]);
            N = length(d);
            t = (0:1/fs: (N-1)/fs) + (s_idx-1)*N/fs;
            dfilt=filtfilt(af,bf,d(:,ch_all(ch_idx)+1));
            plot(t(1:10:end),dfilt(1:10:end),'b','linewidth',1); hold on
        end
        xlim([offset offset+0.4])
        ylim([-1500 1500])
        set(gca,'fontsize',28,'xtick',0:.05:5,'xticklabel',(0:.05:5)-offset,'linewidth',1)
        set(gcf,'position',[0 360 1440 360])
        xlabel('Time (Sec)')
        ylabel('A.U.')
        %title(['BD' num2str(bd_all(ch_idx)) ', Ch' num2str(ch_all(ch_idx))])
        %saveas(gcf,[DirS '/bd_cp30_ds10_' num2str(test) '.tdms_PXI1Slot' num2str(bd_all(ch_idx)) '_ch' num2str(ch_all(ch_idx)) '.png'],'png')
    end
end

%% plot xcorr
j=16;
d_align=wfm_align_bc1f(:,:,j);
if j == 5
    d_align = -d_align; % sensor 5's polarity is reversed
end
Ppick = Ppick_AIC(j); % P pick with AIC algorithm for 1bc bounce
% cut 400 micron, 100 before P and 300 after P
LWin = round(400e-6*fs)+1;
LBefore = round(100e-6*fs);
LAfter = round(300e-6*fs);
s1 = d_align(Ppick - LBefore:round(Ppick + LAfter),1);
[acor,lag] = xcorr(dfilt,s1);