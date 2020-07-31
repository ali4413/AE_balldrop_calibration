% %% write AIC input waveforms
% %load('/Users/chengu/Desktop/AE_test/Berea_SSt_3_31_2016/STALTA_all/STALTA_unfilt_2_50_1.3_1.1.mat')
% load('/Users/chengu/Desktop/AE_test/Berea_SSt_3_31_2016/STALTA_all/STALTA_filt_2_50_1.5_1.1.mat')
% Dir = '/Users/chengu/Desktop/AE_test/Berea_SSt_3_31_2016/Amat/';
% NSeg = 40697728;
% fs=25e6;
% DELTA=1/fs; % 25MHz
% factors = factor(NSeg);
% Twin = prod(factors(1:9));
% Nwin = NSeg/Twin;
% idx = ones(8,1);
% Ntrigs = [length(trigt_all_all{1});length(trigt_all_all{2});length(trigt_all_all{3});length(trigt_all_all{4});...
%     length(trigt_all_all{5});length(trigt_all_all{6});length(trigt_all_all{7});length(trigt_all_all{8})];
% ev_trigt_STALTA_all=nan(sum(Ntrigs),8);
% ev_detrigt_STALTA_all=nan(sum(Ntrigs),8);
% ev_trigt_AIC_all=nan(sum(Ntrigs),8);
% ev_trigt_STALTA=nan(8,2);
% ev_trigt_AIC=nan(8,1);
% ev_wfm=cell(sum(Ntrigs),1);
% ev=0;
% %%
% tic
% while(idx(1)<=Ntrigs(1)||idx(2)<=Ntrigs(2)||idx(3)<=Ntrigs(3)||idx(4)<=Ntrigs(4)|| ...
%         idx(5)<=Ntrigs(5)||idx(6)<=Ntrigs(6)||idx(7)<=Ntrigs(7)||idx(8)<=Ntrigs(8))
%     ev=ev+1;
%     ev_trigt_STALTA=nan(8,2);
%     ev_trigt_AIC=nan(8,1);
%     trigt_tmp=nan(8,2);
%     for ch=0:7
%         if idx(ch+1)<=Ntrigs(ch+1)
%             trigt_tmp(ch+1,:)=trigt_all_all{ch+1}(idx(ch+1),:);
%         end
%     end
%     [trigt_min,ch_idx]=min(trigt_tmp(:,1));
%     trigs_min=round(trigt_min*fs)+1;
%     seg_id = ceil(trigs_min/NSeg);
%     ev_trigt_STALTA(ch_idx,:)=trigt_tmp(ch_idx,:);
%     idx(ch_idx)=idx(ch_idx)+1;
%     idx_same_ev=find(trigt_tmp(:,1)-trigt_min>0 & trigt_tmp(:,1)-trigt_min<20e-6);
%     ev_trigt_STALTA(idx_same_ev,:)=trigt_tmp(idx_same_ev,:);
%     idx(idx_same_ev)=idx(idx_same_ev)+1;
%     %     %%
%     %     N=ev_Es-ev_Bs+1;
%     %     wfm_temp=zeros(N,8);
%     %     for ch=0:7
%     %         testfile = ['AE_ch' num2str(ch+1) '_' num2str(seg_id) '.mat'];
%     %         load([Dir testfile]);
%     %         Rw=d(ev_Bs:ev_Es);
%     %         wfm_temp(:,ch+1)=Rw;
%     %         AIC=zeros(N,1);
%     %         for tw=2:N-1
%     %             %if var(Rw(1:tw))==0 || var(Rw(tw+1:end))==0
%     %             %   AIC(tw)=0;
%     %             %else
%     %             AIC(tw)=tw*log(var(Rw(1:tw)))+(N-tw-1)*log(var(Rw(tw+1:end)));
%     %             %end
%     %         end
%     %         AIC(1)=AIC(2);
%     %         AIC(end)=AIC(end-1);
%     % %         figure;
%     % %         subplot(2,1,1)
%     % %         plot(Rw)
%     % %         subplot(2,1,2)
%     % %         plot(AIC)
%     %         [tmp,Ptime_pick_AIC]=min(AIC);
%     %         ev_trigt_AIC(ch+1)=Ptime_pick_AIC;
%     %     end
%     %     %%
%     ev_trigt_STALTA_all(ev,:)=ev_trigt_STALTA(:,1)';
%     ev_detrigt_STALTA_all(ev,:)=ev_trigt_STALTA(:,2)';
%     %ev_trigt_AIC_all(ev,:)=ev_trigt_AIC;
%     %ev_wfm{ev}=wfm_temp;
% end
% toc
% % %%
% save STALTA_ev_8ch_filt_10k_500k.mat
% %% 8 channel record for detected events
% tic
% load STALTA_ev_8ch.mat
%num=ev;
seg_id0=0;
d8=zeros(NSeg,8);
for ev=620970
    ev_B=max([1 min(min(ev_trigt_STALTA_all(ev,:)))-50e-6]);
    ev_E=max(max(ev_detrigt_STALTA_all(ev,:)))+50e-6;
    trigs_B=round(ev_B*fs)+1;
    seg_id = ceil(trigs_B/NSeg);
    ev_Bs=round(ev_B*fs)-NSeg*(seg_id-1);
    ev_Es=round(ev_E*fs)-NSeg*(seg_id-1);
    N=ev_Es-ev_Bs+1;
    if seg_id>seg_id0
        for ch=0:7
            testfile = ['/Amat/AE_ch' num2str(ch+1) '_' num2str(seg_id) '.mat'];
            load([Dir testfile]);
            d8(:,ch+1)=d';
        end
    end
    Rw=d8(ev_Bs:ev_Es,:);
    wfm_temp=Rw;
%     AIC=zeros(N,8);
%     for tw=2:N-1
%         %if var(Rw(1:tw))==0 || var(Rw(tw+1:end))==0
%         %   AIC(tw)=0;
%         %else
%         AIC(tw,:)=tw*log(var(Rw(1:tw,:)))+(N-tw-1)*log(var(Rw(tw+1:end,:)));
%         %end
%     end
%     AIC(1,:)=AIC(2,:);
%     AIC(end,:)=AIC(end-1,:);
%     %         figure;
%     %         subplot(2,1,1)
%     %         plot(Rw)
%     %         subplot(2,1,2)
%     %         plot(AIC)
%     [tmp,Ptime_pick_AIC]=min(AIC);
%     seg_id0=seg_id;
%     ev_trigt_AIC_all(ev,:)=Ptime_pick_AIC;
%     ev_wfm{ev}=wfm_temp;
%     if mod(ev,100)==0
%         display(ev)
%         save AIC_test_wfm_10k_500k.mat
%     end
end
% toc
% %%
% for ev=1:num
%     h=figure;
%     for i=0:7
%         subplot(16,1,i*2+1:i*2+2)
%         t=0:1/25e6:(length(ev_wfm{ev}(:,i+1))-1)/25e6;
%         plot(t,ev_wfm{ev}(:,i+1),'linewidth',1)
%         xlim([t(1) t(end)])
%         ylabel(['ch ' num2str(i)],'fontsize',14,'fontname','arial')
%         if i==0
%             title(['EV ' num2str(ev)])
%         end
%         if i==7
%             xlabel('Time (sec)','fontsize',14,'fontname','arial')
%         else
%             set(gca,'xticklabel',[])
%         end
%         set(gca,'fontsize',14,'fontname','arial')
%         set(gcf,'Position',[30 30 500 1000])
%     end
%     saveas(h,['/nobackup1b/users/guchch/AE_test/Berea_SSt_3_31_2016/AIC_all/fig_wfm/wfm_ev_' num2str(ev) '.fig'],'fig')
%     saveas(h,['/nobackup1b/users/guchch/AE_test/Berea_SSt_3_31_2016/AIC_all/fig_wfm/wfm_ev_' num2str(ev) '.eps'],'eps')
%     close(h)
% end