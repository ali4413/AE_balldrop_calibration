for ch_idx = 1
    %%
    r_vec = r_vec_all(ch_idx,:);
    Ppick = Ppick_AIC(ch_idx); % P pick with AIC algorithm for 1bc bounce
    d_align=wfm_align_bc1f(:,:,ch_idx);
    if ch_idx ==5
        d_align=-d_align;
    end
    d_obs = reshape(d_align(Ppick - LBefore:round(Ppick + LAfter),:),LWin*3,1);
    d_obs_norm = d_obs/norm(d_obs);
    burnIn = 2e4;
    %matfile = ['MCMC_4para_sensor_' num2str(ch_idx) '_PpickAIC1bc.mat'];
    matfile = ['results_mat/MCMC_4para_scale_sigma_sensor_' num2str(ch_idx) '_nlev' num2str(nlev_all(ch_idx)) '_PpickAIC1bc.mat'];
    load(matfile);
    wfm_post=zeros(15003,round((Nm-burnIn)/10)+1);
    for i=burnIn:10:Nm
        [syn1, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0(i),omega_s(i), varepsilon(i));
        [syn2, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0(i)*a(i),omega_s(i), varepsilon(i));
        [syn3, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0(i)*a(i)^2,omega_s(i), varepsilon(i));
        d_align_syn  = [syn1 syn2 syn3];
        Tp_sample = round(sqrt(sum(r_vec.^2))/Vp*fs/1e3); %samples before P arrival
        d_syn_tmp = [zeros(LBefore,3); d_align_syn(Tp_sample:Tp_sample+LAfter,:)];
        d_syn = reshape(filtfilt(af,bf, d_syn_tmp),LWin*3,1);
        d_syn_norm = d_syn/norm(d_syn);
        wfm_post(:,(i-burnIn)/10+1) = d_syn_norm;
    end
    %%
    h1=figure;
    h1name = ['wfm_' num2str(ch_idx) '_PpickAIC1bc.png'];
    tmicron = 1e6*(0:dt:5e3*dt);
    subplot(3,1,1)
    plot(tmicron,d_obs_norm(1:LWin),'r','linewidth',1); hold on
    plot(tmicron,wfm_post(1:LWin,:),'color',[.5, .5, .5],'linewidth',1);
    plot(tmicron,mean(wfm_post(1:LWin,:),2),'b','linewidth',1);
    ylim([-.08 .08])
    set(gca,'fontsize',14,'xticklabel',[])
    title('1st Bounce')
    legend('Observed waveforms','Mean posterior predicted waveforms')
    subplot(3,1,2)
    plot(tmicron,d_obs_norm(LWin+1:LWin*2),'r','linewidth',1); hold on
    plot(tmicron,wfm_post(LWin+1:LWin*2,:),'color',[.5, .5, .5],'linewidth',1)
    plot(tmicron,mean(wfm_post(LWin+1:LWin*2,:),2),'b','linewidth',1);
    ylim([-.08 .08])
    set(gca,'fontsize',14,'xticklabel',[])
    title('2nd Bounce')
    subplot(3,1,3)
    plot(tmicron,d_obs_norm(LWin*2+1:LWin*3),'r','linewidth',1); hold on
    plot(tmicron,wfm_post(LWin*2+1:LWin*3,:),'color',[.5, .5, .5],'linewidth',1)
    plot(tmicron,mean(wfm_post(LWin*2+1:LWin*3,:),2),'b','linewidth',1);
    ylim([-.08 .08])
    set(gca,'fontsize',14)
    title('3rd Bounce')
    xlabel('Time (\muS)')
    set(h1,'position',[0 0 560 1024])
end 
