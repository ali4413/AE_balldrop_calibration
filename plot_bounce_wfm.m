DirM = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc/mat_bd/bd_cp30_ds10';
DirS = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc/wfm_plot_bounce';
DirD = '/Users/chengu/Desktop/AE_calibration/AE_mat_balldrop_16ch_@pc';
%load([DirD 'wfm_bounce.mat'])
% fs = 12.5e6;
% bd_all = [repmat(6,1,8) repmat(8,1,8)];
% ch_all = repmat(0:7,1,2);
%t_detect = cell(16,1);
% %%
% figure;
% for i=1:16
%     plot(t_detect{i},'o'); hold on
% end
% [t_bounce_rough,y]=ginput(7);
% %%
% t_bounce_rough = [3.3 6.9 10.1 13.4 16.7 19.5 22.1];
% Nbounce = length(t_bounce_rough);
N = 8854272; % length of one data segment
% t_bounce = zeros(16,21);
% B_bounce = zeros(7,3);
% wfm_bounce = zeros(20001,16,21);
for test = 3
    for bc = 1:Nbounce
        %%
        for ch_idx = 1:16
            tmp = t_detect{ch_idx};
            t_bounce_tmp = uniquetol(tmp(tmp>t_bounce_rough(bc)-.5&tmp<t_bounce_rough(bc)+.5),1e-5);
            for i = 1:length(t_bounce_tmp)
                t_bounce(ch_idx,bc+(i-1)*7)=t_bounce_tmp(i);
            end
        end
        %%
        for i = 1:3
            tmp = t_bounce(t_bounce(:,bc+(i-1)*7)~=0,bc+(i-1)*7);
            if tmp~=0
                if i==2 && bc==2 || bc==1 && i==3
                    B_bounce(bc,i) = median(tmp); % there is an detected time outlier
                else
                    B_bounce(bc,i) = min(tmp);
                end
                
                s_idx = ceil((round(B_bounce(bc,i)*fs)+1)/N);
                Bs = round(B_bounce(bc,i)*fs)+1 - N*(s_idx-1);
                %if i==1
                filename = ['bd_cp30_ds10_' num2str(test) '.tdms_PXI1Slot' num2str(6) '_' num2str(s_idx) '.mat'];
                load([DirM '/' filename],'d');
                d6=d;
                %end
                wfm_bounce(:,1:8,bc+7*(i-1)) = d6(Bs-1e4:Bs+1e4,:);
                %if i==1
                filename = ['bd_cp30_ds10_' num2str(test) '.tdms_PXI1Slot' num2str(8) '_' num2str(s_idx) '.mat'];
                load([DirM '/' filename],'d');
                d8=d;
                %end
                wfm_bounce(:,9:16,bc+7*(i-1)) = d8(Bs-1e4:Bs+1e4,:);
            end
        end
    end
end

%% plot
amax=max(abs(wfm_bounce(:)));
for i = 1:3
    for bc = 1
        filename=['wfm_bounce_' num2str(bc) '_' num2str(i)];
        h=figure; hold on;
        for ch_idx = 1:16
            if ch_idx<=8
                hold on
                subplot('position',[0.12+l_d(ch_idx)/200 -.05+Axial_d(ch_idx)/80 0.8/10 1.0/5]); hold on
                plot(0:1/fs*1e6:2e4/fs*1e6,wfm_bounce(:,ch_idx,bc+7*(i-1))/amax,'k','linewidth',1)
                set(gca,'fontsize',18)
                title(['BD 6, Ch ' num2str(ch_all(ch_idx)) '.'])
                ylim([-1 1])
                xlim([700 1e3])
                axis off
            else
                hold on
                subplot('position',[0.12+l(ch_idx-8)/200 -.05+Axial(ch_idx-8)/80 0.8/10 1.0/5]); hold on
                plot(0:1/fs*1e6:2e4/fs*1e6,wfm_bounce(:,ch_idx,bc+7*(i-1))/amax,'r','linewidth',1)
                set(gca,'fontsize',18)
                title(['BD 8, Ch ' num2str(ch_all(ch_idx))])
                ylim([-1 1])
                xlim([700 1e3])
                axis off
            end
        end
        set(gcf,'position',[0 0 1440 960])
        %saveas(h,[DirS '/' filename '.fig'],'fig')
        %saveas(h,[DirS '/' filename '.png'],'png')
    end
end