DirM = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc/mat_bd/bd_cp30_ds10';
DirS = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc/wfm_plot_continous';
DirD = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc';
fs = 12.5e6;
bd_all = [repmat(6,1,8) repmat(8,1,8)];
ch_all = repmat(0:7,1,2);
%t_detect = cell(16,1);
for test = 3
    load([DirD '/bd_cp30_ds10_' num2str(test) '.tdms_stalta.mat'])
    for ch_idx = [1:7 9:15]
        h=figure;
        for s_idx = 1 : 32
            filename = ['bd_cp30_ds10_' num2str(test) '.tdms_PXI1Slot' num2str(bd_all(ch_idx)) '_' num2str(s_idx) '.mat'];
            load([DirM '/' filename]);
            N = length(d);
            t = (0:1/fs: (N-1)/fs) + (s_idx-1)*N/fs;
            plot(t(1:10:end),d(1:10:end,ch_all(ch_idx)+1),'b','linewidth',1); hold on
            if ch_idx <= 8
                trigs_tmp = trigs_all{s_idx,ch_all(ch_idx)+1};
            else
                trigs_tmp = trigs_all{s_idx+32,ch_all(ch_idx)+1};
            end
            if ~isempty(trigs_tmp)
                lt = length(trigs_tmp);
                for n = 1:lt
                    [dmax, smax] = max(d((trigs_tmp(n)-3e4):(trigs_tmp(n)+3e4),ch_all(ch_idx)+1));
                    plot(t(trigs_tmp(n)-3e4-1+smax),dmax+100,'v','markersize',12,'markerfacecolor','r','markeredgecolor','r')
                    t_detect{ch_idx} = [t_detect{ch_idx} t(trigs_tmp(n)-3e4-1+smax)];
                end
            end
        end
        xlim([0 N*32/fs])
        set(gca,'fontsize',28,'xtick',0:2:25,'linewidth',1)
        set(gcf,'position',[0 360 1440 360])
        xlabel('Time (Sec)')
        ylabel('A.U.')
        title(['BD' num2str(bd_all(ch_idx)) ', Ch' num2str(ch_all(ch_idx))])
        %saveas(gcf,[DirS '/bd_cp30_ds10_' num2str(test) '.tdms_PXI1Slot' num2str(bd_all(ch_idx)) '_ch' num2str(ch_all(ch_idx)) '.png'],'png')
    end
end