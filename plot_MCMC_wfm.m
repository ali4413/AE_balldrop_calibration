mu=zeros(4,16,9);
sigma=zeros(4,16,9);
for ch_idx = [1:12 14:16]
    for nlev=.1:.1:.9
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
        matfile = ['results_mat/MCMC_4para_scale_sigma_sensor_' num2str(ch_idx) '_nlev' num2str(nlev) '_PpickAIC1bc.mat'];
        load(matfile);
        %
        [syn1, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,mean(v0(burnIn+1:end)),mean(omega_s(burnIn+1:end)), mean(varepsilon(burnIn+1:end)));
        [syn2, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,mean(v0(burnIn+1:end))*mean(a(burnIn+1:end)),mean(omega_s(burnIn+1:end)), mean(varepsilon(burnIn+1:end)));
        [syn3, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2,mean(omega_s(burnIn+1:end)), mean(varepsilon(burnIn+1:end)));
        d_align_syn  = [syn1 syn2 syn3];
        Tp_sample = round(sqrt(sum(r_vec.^2))/Vp*fs/1e3); %samples before P arrival
        d_syn_tmp = [zeros(LBefore,3); d_align_syn(Tp_sample:Tp_sample+LAfter,:)];
        d_syn = reshape(filtfilt(af,bf, d_syn_tmp),LWin*3,1);
        d_syn_norm = d_syn/norm(d_syn);
        h1=figure;
        h1name = ['results_fig/wfm_' num2str(ch_idx) '_nlev' num2str(nlev) '_PpickAIC1bc.eps'];
        tmicron = 1e6*(0:dt:5e3*dt);
        subplot(3,1,1)
        plot(tmicron,d_obs_norm(1:LWin),'r','linewidth',1); hold on
        plot(tmicron,d_syn_norm(1:LWin),'b','linewidth',1);
        text(10,.01,'1st Bounce','fontsize',14)
        ylim([-.08 .08])
        set(gca,'fontsize',14,'xticklabel',[])
        title('1st Bounce')
        legend('Observed waveforms','Mean posterior predicted waveforms')
        subplot(3,1,2)
        plot(tmicron,2*mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/2/g+d_obs_norm(LWin+1:LWin*2),'r','linewidth',1); hold on
        plot(tmicron,2*mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/2/g+d_syn_norm(LWin+1:LWin*2),'b','linewidth',1)
        ylim([-.08 .08])
        set(gca,'fontsize',14)
        text(10,.09,'2nd Bounce','fontsize',14)
        plot(tmicron,2*(mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/2/g+mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2/2/g)+d_obs_norm(LWin*2+1:LWin*3),'r','linewidth',1); hold on
        plot(tmicron,2*(mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/2/g+mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2/2/g)+d_syn_norm(LWin*2+1:LWin*3),'b','linewidth',1)
        ylim([-.08 .08])
        set(gca,'fontsize',14)
        text(10,.14,'3rd Bounce','fontsize',14)
        xlabel('Time (\muS)')
        axis off
        set(h1,'position',[0 0 560 1024])
        saveas(h1,h1name,'epsc')
        close(h1)
        

        mu(:,ch_idx,round(nlev/.1))=[mean(v0(burnIn+1:end));mean(a(burnIn+1:end));mean(omega_s(burnIn+1:end));mean(varepsilon(burnIn+1:end))];
        sigma(:,ch_idx,round(nlev/.1))=[std(v0(burnIn+1:end));std(a(burnIn+1:end));std(omega_s(burnIn+1:end));std(varepsilon(burnIn+1:end))];
        h2=figure;
        h2name = ['results_fig/MCMC_chain_2_' num2str(ch_idx) '_nlev' num2str(nlev) '_PpickAIC1bc.png'];
        subplot(4,1,1); plot(v0);
        set(gca,'fontsize',14,'xticklabel',[])
        ylim([mean(v0(burnIn+1:end))-3*std(v0(burnIn+1:end)) mean(v0(burnIn+1:end))+3*std(v0(burnIn+1:end))])
        title('Initial velocity')
        ylabel('v_0 (m/s)')
        subplot(4,1,2); plot(a);
        set(gca,'fontsize',14,'xticklabel',[])
        ylim([mean(a(burnIn+1:end))-3*std(a(burnIn+1:end)) mean(a(burnIn+1:end))+3*std(a(burnIn+1:end))])
        title('Rebound coefficient')
        ylabel('a')
        subplot(4,1,3); plot(omega_s/1e3);
        set(gca,'fontsize',14,'xticklabel',[])
        ylim([mean(omega_s(burnIn+1:end))-3*std(omega_s(burnIn+1:end)) mean(omega_s(burnIn+1:end))+3*std(omega_s(burnIn+1:end))]/1e3)
        title('Resonance frequency')
        ylabel('\omega_s (kHz)')
        subplot(4,1,4); plot(varepsilon/1e3);
        ylim([mean(varepsilon(burnIn+1:end))-3*std(varepsilon(burnIn+1:end)) mean(varepsilon(burnIn+1:end))+3*std(varepsilon(burnIn+1:end))]/1e3)
        set(gca,'fontsize',14)
        title('Damping coefficient')
        ylabel('\epsilon (kHz)')
        xlabel('MCMC iterations')
        saveas(h2,h2name,'png')
        close(h2)
        %%
        h3=figure;
        h3name = ['results_fig/MCMC_matrix2_' num2str(ch_idx) '_nlev' num2str(nlev) '_PpickAIC1bc.png'];
        [H,AX,BigAx,P,PAx]=plotmatrix([v0(burnIn+1:end) a(burnIn+1:end) omega_s(burnIn+1:end)/1e3 varepsilon(burnIn+1:end)/1e3]);
        labeltext={'v_0 (m/s)','a', '\omega_s (kHz)','\epsilon (kHz)'};
        for i=1:4
            for j=1:4
                set(AX(i,j),'fontsize',14)
            end
            ylabel(AX(i),labeltext{i},'fontsize',14)
            xlabel(AX(i*4),labeltext{i},'fontsize',14)
        end
        set(gca,'fontsize',14,'fontname','Arial')
        saveas(h3,h3name,'png')
        close(h3)
        %% trajectory
        h4 = figure;
        h4name = ['results_fig/Trajectory2_' num2str(ch_idx) '_nlev' num2str(nlev) '_PpickAIC1bc.png'];
        
        % observation
        tt=(0:1e-3:.45)';
        height = zeros(length(tt),9e3);
        for k = burnIn+1:Nm
            t1 = v0(k)/g;
            t2 = v0(k)*a(k)/g*2;
            t3 = v0(k)*a(k)^2/g*2;
            h_tmp = [v0(k)^2/2/g-1/2*g*tt(tt<t1).^2; ...
                v0(k)*a(k)*(tt(tt>t1&tt<t1+t2)-t1)-1/2*g*(tt(tt>t1&tt<t1+t2)-t1).^2; ...
                v0(k)*a(k)^2*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2)-1/2*g*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2).^2;...
                v0(k)*a(k)^3*(tt(tt>t1+t2+t3)-t1-t2-t3)-1/2*g*(tt(tt>t1+t2+t3)-t1-t2-t3).^2];
            height(:,k-burnIn)=h_tmp;
        end
        
        t1 = mean(v0(burnIn+1:end))/g;
        t2 = mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/g*2;
        t3 = mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2/g*2;
        mean_height = [mean(v0(burnIn+1:end))^2/2/g-1/2*g*tt(tt<t1).^2; ...
            mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))*(tt(tt>t1&tt<t1+t2)-t1)-1/2*g*(tt(tt>t1&tt<t1+t2)-t1).^2; ...
            mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2)-1/2*g*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2).^2;...
            mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^3*(tt(tt>t1+t2+t3)-t1-t2-t3)-1/2*g*(tt(tt>t1+t2+t3)-t1-t2-t3).^2];
        
        plot(height(:,1:100:end),tt,'color',[.5 .5 .5]); hold on
        plot(mean_height,tt,'r','linewidth',1)
        
        ylabel('Bouncing time interval (sec)')
        xlabel('Height (m)')
        title('Posterior Predicted Trajectory')
        xlim([0 .12])
        ylim([tt(1) tt(end)])
        set(gca,'fontsize',18,'xdir','reverse')
        saveas(h4,h4name,'png')
        close(h4)
    end
end
