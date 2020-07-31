res_table = zeros(16,11);
Nm=1e6;
for ch_idx = [1:12 14:16]
    %for ch_idx = 16
    %%
    %res_table(ch_idx,:)=[]
    r_vec = r_vec_all(ch_idx,:);
    Ppick = Ppick_AIC(ch_idx); % P prick with AIC algorithm for 1bc bounce
    d_align=wfm_align_bc1f(:,:,ch_idx);
    if ch_idx ==5
        d_align=-d_align;
    end
    d_obs = reshape(d_align(Ppick - LBefore:round(Ppick + LAfter),:),LWin*3,1);
    d_obs_norm = d_obs/norm(d_obs);
    alpha=sum(d_obs_norm(1:LBefore).^2)/LBefore*3*LWin;
    burnIn = 6e5;
    matfile = ['results_mat_9/MCMC_4para_scale_sigma_sensor_' num2str(ch_idx) '_k_1000000_PpickAIC1bc_9.mat'];
    wfmfile = ['results_mat_9/wfm_' num2str(ch_idx) '.mat'];
    load(matfile);
    load(wfmfile);
    %     res_table(ch_idx,:)=[alpha mean(sdt(burnIn+1:end)) std(sdt(burnIn+1:end)) ...
    %         mean(v0(burnIn+1:end)) std(v0(burnIn+1:end))...
    %         mean(a(burnIn+1:end)) std(a(burnIn+1:end))...
    %         mean(omega_s(burnIn+1:end))/1e3 std(omega_s(burnIn+1:end))/1e3...
    %         mean(varepsilon(burnIn+1:end))/1e3 std(varepsilon(burnIn+1:end))/1e3];
    %%
    wfm_mean=mean(wfm_test,2);
    wfm_std=std(wfm_test,[],2);
    h1 = figure;
    h1name = ['results_fig_9/wfm_' num2str(ch_idx) '_notext_blue.eps'];
    tmicron = 1e6*(0:dt:5e3*dt);
    X1=[tmicron,(fliplr(tmicron))];                %#create continuous x value array for plotting
    Y12=[(wfm_mean(1:LWin) - 2*wfm_std(1:LWin))',fliplr((wfm_mean(1:LWin) + 2*wfm_std(1:LWin))')];           %#create y values for out and then back
    h2bc=fill(X1,Y12,[0.6353    0.8118    0.9961]); hold on
    set(h2bc,'EdgeColor','none')
    Y11=[(wfm_mean(1:LWin) - wfm_std(1:LWin))',fliplr((wfm_mean(1:LWin) + wfm_std(1:LWin))')];              %#create y values for out and then back
    h1bc=fill(X1,Y11,[0.0118    0.2627    0.8745]);
    set(h1bc,'EdgeColor','none')
    plot(tmicron,d_obs_norm(1:LWin),'r','linewidth',1);
    plot(tmicron,wfm_mean(1:LWin),'k','linewidth',1);
    %text(10,.01,'1st Bounce','fontsize',14)
    
    offset1=mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/g;
    X2=[tmicron,(fliplr(tmicron))];                %#create continuous x value array for plotting
    Y22=[(offset1+wfm_mean(LWin+1:LWin*2) - 2*wfm_std(LWin+1:LWin*2))',fliplr(offset1+(wfm_mean(LWin+1:LWin*2) + 2*wfm_std(LWin+1:LWin*2))')];           %#create y values for out and then back
    h2bc=fill(X1,Y22,[0.6353    0.8118    0.9961]); hold on
    set(h2bc,'EdgeColor','none')
    Y21=[(offset1+wfm_mean(LWin+1:LWin*2) - wfm_std(LWin+1:LWin*2))',fliplr(offset1+(wfm_mean(LWin+1:LWin*2) + wfm_std(LWin+1:LWin*2))')];              %#create y values for out and then back
    h1bc=fill(X1,Y21,[0.0118    0.2627    0.8745]);
    set(h1bc,'EdgeColor','none')
    plot(tmicron,offset1+d_obs_norm(LWin+1:LWin*2),'r','linewidth',1); hold on
    plot(tmicron,offset1+wfm_mean(LWin+1:LWin*2),'k','linewidth',1)
    set(gca,'fontsize',14)
    %text(10,offset1+.01,'2nd Bounce','fontsize',14)
    offset2=(mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/g+mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2/g);
    X3=[tmicron,(fliplr(tmicron))];                %#create continuous x value array for plotting
    Y32=[(offset2+wfm_mean(LWin*2+1:LWin*3) - 2*wfm_std(LWin*2+1:LWin*3))',fliplr((offset2+wfm_mean(LWin*2+1:LWin*3) + 2*wfm_std(LWin*2+1:LWin*3))')];           %#create y values for out and then back
    h2bc=fill(X1,Y32,[0.6353    0.8118    0.9961]); hold on
    set(h2bc,'EdgeColor','none')
    Y31=[(offset2+wfm_mean(LWin*2+1:LWin*3) - wfm_std(LWin*2+1:LWin*3))',fliplr((offset2+wfm_mean(LWin*2+1:LWin*3) + wfm_std(LWin*2+1:LWin*3))')];              %#create y values for out and then back
    h1bc=fill(X1,Y31,[0.0118    0.2627    0.8745]);
    set(h1bc,'EdgeColor','none')
    plot(tmicron,offset2+d_obs_norm(LWin*2+1:LWin*3),'r','linewidth',1); hold on
    plot(tmicron,offset2+wfm_mean(LWin*2+1:LWin*3),'k','linewidth',1)
    set(gca,'fontsize',14)
    %text(10,offset2+.01,'3rd Bounce','fontsize',14)
    xlabel('Time (\muS)')
    ylim([-mean(v0(burnIn+1:end))/g/2 .225-mean(v0(burnIn+1:end))/g/2])
    axis off
    saveas(h1,h1name,'epsc')
    %     close(h1)
    %     %%
    %     h2=figure;
    %     h2name = ['results_fig_6/MCMC_chain_' num2str(ch_idx) '_PpickAIC1bc_2.png'];
    %     subplot(8,1,1:2); plot(v0);
    %     set(gca,'fontsize',14,'xticklabel',[])
    %     ylim([mean(v0(burnIn+1:end))-3*std(v0(burnIn+1:end)) mean(v0(burnIn+1:end))+3*std(v0(burnIn+1:end))])
    %     %title('Initial velocity')
    %     ylabel('v_0 (m/s)')
    %     subplot(8,1,3:4); plot(a);
    %     set(gca,'fontsize',14,'xticklabel',[])
    %     ylim([mean(a(burnIn+1:end))-3*std(a(burnIn+1:end)) mean(a(burnIn+1:end))+3*std(a(burnIn+1:end))])
    %     %title('Rebound coefficient')
    %     ylabel('a')
    %     subplot(8,1,5:6); plot(omega_s/1e3);
    %     set(gca,'fontsize',14,'xticklabel',[])
    %     ylim([mean(omega_s(burnIn+1:end))-3*std(omega_s(burnIn+1:end)) mean(omega_s(burnIn+1:end))+3*std(omega_s(burnIn+1:end))]/1e3)
    %     %title('Resonance frequency')
    %     ylabel('\omega_s (kHz)')
    %     subplot(8,1,7:8); plot(varepsilon/1e3);
    %     ylim([mean(varepsilon(burnIn+1:end))-3*std(varepsilon(burnIn+1:end)) mean(varepsilon(burnIn+1:end))+3*std(varepsilon(burnIn+1:end))]/1e3)
    %     set(gca,'fontsize',14)
    %     %title('Damping coefficient')
    %     ylabel('\epsilon (kHz)')
    %     xlabel('MCMC iterations')
    %     saveas(h2,h2name,'png')
    %     %close(h2)
    %     %%
    %     h3=figure;
    %     h3name = ['results_fig_6/MCMC_matrix_' num2str(ch_idx) '_PpickAIC1bc_2.png'];
    %     [H,AX,BigAx,P,PAx]=plotmatrix([v0(burnIn+1:end) a(burnIn+1:end) omega_s(burnIn+1:end)/1e3 varepsilon(burnIn+1:end)/1e3]);
    %     labeltext={'v_0 (m/s)','a', '\omega_s (kHz)','\epsilon (kHz)'};
    %     for i=1:4
    %         for j=1:4
    %             set(AX(i,j),'fontsize',14)
    %         end
    %         ylabel(AX(i),labeltext{i},'fontsize',14)
    %         xlabel(AX(i*4),labeltext{i},'fontsize',14)
    %     end
    %     set(gca,'fontsize',14,'fontname','Arial')
    %     saveas(h3,h3name,'png')
    %     %close(h3)
    %     %% trajectory
    %     h4 = figure;
    %     h4name = ['results_fig_6/Trajectory_' num2str(ch_idx) '_PpickAIC1bc_2.png'];
    %
    %     % observation
    %     tt=(0:1e-3:.45)';
    %     height = zeros(length(tt),9e3);
    %     for k = burnIn+1:Nm
    %         t1 = v0(k)/g;
    %         t2 = v0(k)*a(k)/g*2;
    %         t3 = v0(k)*a(k)^2/g*2;
    %         h_tmp = [v0(k)^2/2/g-1/2*g*tt(tt<t1).^2; ...
    %             v0(k)*a(k)*(tt(tt>t1&tt<t1+t2)-t1)-1/2*g*(tt(tt>t1&tt<t1+t2)-t1).^2; ...
    %             v0(k)*a(k)^2*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2)-1/2*g*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2).^2;...
    %             v0(k)*a(k)^3*(tt(tt>t1+t2+t3)-t1-t2-t3)-1/2*g*(tt(tt>t1+t2+t3)-t1-t2-t3).^2];
    %         height(:,k-burnIn)=h_tmp;
    %     end
    %
    %     t1 = mean(v0(burnIn+1:end))/g;
    %     t2 = mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/g*2;
    %     t3 = mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2/g*2;
    %     mean_height = [mean(v0(burnIn+1:end))^2/2/g-1/2*g*tt(tt<t1).^2; ...
    %         mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))*(tt(tt>t1&tt<t1+t2)-t1)-1/2*g*(tt(tt>t1&tt<t1+t2)-t1).^2; ...
    %         mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2)-1/2*g*(tt(tt>t1+t2&tt<t1+t2+t3)-t1-t2).^2;...
    %         mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^3*(tt(tt>t1+t2+t3)-t1-t2-t3)-1/2*g*(tt(tt>t1+t2+t3)-t1-t2-t3).^2];
    %
    %     plot(height(:,1:100:end),tt,'color',[.5 .5 .5]); hold on
    %     plot(mean_height,tt,'r','linewidth',1)
    %
    %     ylabel('Bouncing time interval (sec)')
    %     xlabel('Height (m)')
    %     title('Posterior Predicted Trajectory')
    %     xlim([0 .12])
    %     ylim([tt(1) tt(end)])
    %     set(gca,'fontsize',18,'xdir','reverse')
    %     saveas(h4,h4name,'png')
    %     %close(h4)
    
end
