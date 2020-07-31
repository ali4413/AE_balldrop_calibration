%%
ev_trigt_AIC_all=zeros(2,8);
for G=1:2
    for ev=1
        Rw=ev_wfm{G,ev};
        N=length(Rw);
        AIC=zeros(N,8);
        for tw=2:N-1
            %if var(Rw(1:tw))==0 || var(Rw(tw+1:end))==0
            %   AIC(tw)=0;
            %else
            AIC(tw,:)=tw*log(var(Rw(1:tw,:)))+(N-tw-1)*log(var(Rw(tw+1:end,:)));
            %end
        end
        AIC(1,:)=AIC(2,:);
        AIC(end,:)=AIC(end-1,:);
        AIC_range=100:1.25e3;
        [tmp,Ptime_pick_AIC]=min(AIC(AIC_range,:));
        Ptime_pick_AIC=Ptime_pick_AIC+AIC_range(1)-1;
        %seg_id0=seg_id;
        ev_trigt_AIC_all(G,:)=Ptime_pick_AIC;
        %%
        h=figure;
        t=0:1/fs:1/fs*(N-1);
        for ch=0:7
            subplot(16,2,ch*4+1:2:ch*4+3)
            plot(t*1e6,Rw(:,ch+1),'linewidth',1); hold on
            plot([t(Ptime_pick_AIC(ch+1)) t(Ptime_pick_AIC(ch+1))]*1e6,[-max(abs(Rw(:,ch+1))) max(abs(Rw(:,ch+1)))],'r','linewidth',1)
            xlim([t(1) t(end)]*1e6)
            ylabel(['Ch' num2str(ch)],'fontsize',12)
            if ch==7
                set(gca,'fontsize',10)
                xlabel('Time (\muSec)')
            else
                set(gca,'fontsize',10,'xticklabel',[])
            end
            if ch==0
                title(['PXI1Slot' num2str(board(G))]);
            end
            subplot(16,2,ch*4+2:2:ch*4+4)
            plot(t*1e6,AIC(:,ch+1),'linewidth',1); hold on
            plot([t(Ptime_pick_AIC(ch+1)) t(Ptime_pick_AIC(ch+1))]*1e6,[min(AIC(AIC_range,ch+1)) max(AIC(AIC_range,ch+1))],'r','linewidth',1)
            ylim([min(AIC(:,ch+1)) max(AIC(:,ch+1))])
            xlim([t(1) t(end)]*1e6)
            if ch==7
                set(gca,'fontsize',10)
                xlabel('Time (\muSec)')
            else
                set(gca,'fontsize',10,'xticklabel',[])
            end
            if ch==0
                title('AIC')
            end
            set(gcf,'Position',[0 0 1000 800])
        end
        saveas(h,['/net/quake/archive/acoustic_emission/AE_mat_balldrop_16ch_@pc/AIC_plot/AIC_wfm_' num2str(G) '_' num2str(ev) '.fig'],'fig')
        saveas(h,['/net/quake/archive/acoustic_emission/AE_mat_balldrop_16ch_@pc/AIC_plot/AIC_wfm_' num2str(G) '_' num2str(ev) '.eps'],'epsc')
        %close(h)
    end
end