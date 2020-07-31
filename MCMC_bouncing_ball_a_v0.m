%%
wfm_bc1 = wfm_bounce(:,:,1:7:15);
t_half_bc1 = B_bounce(1,:);
%%
dt = 1/fs;
fn=1/(2*dt);
f1=5e3; f2=500e3;
[a,b]=butter(2,[f1/fn  f2/fn]);
wfm_bc1f=filtfilt(a,b,wfm_bc1);

wfm_crop_bc1f=permute(wfm_bc1f(5e3:15e3,:,:),[1 3 2]);
wfm_align_bc1f=zeros(1e4+1,3,16);
lagDiff = zeros(16,3);
for ch_idx=1:16
    %figure;
    %plot(wfm_crop_bc1f(:,:,ch_idx))
    [cr, lgs]=xcorr(wfm_crop_bc1f(:,:,ch_idx),'coeff');
    for i=1:3
        [~,I] = max(abs(cr(:,i)));
        lagDiff(ch_idx,i) = lgs(I);
        wfm_align_bc1f(:,i,ch_idx)=ddelay(wfm_crop_bc1f(:,i,ch_idx),lagDiff(ch_idx,i));
    end
    %figure;
    %plot(wfm_align_bc1f(:,:,ch_idx)./repmat(max(abs(wfm_align_bc1f(:,:,ch_idx))),1e4+1,1))
end


[X, outliers_idx] = outliers(lagDiff(:,3), 1);
DT=[t_half_bc1(2)-t_half_bc1(1) t_half_bc1(3)-t_half_bc1(2)] ...
    + [lagDiff(~isnan(X),2) lagDiff(~isnan(X),3)-lagDiff(~isnan(X),2)]/fs;

%% MCMC for v0 and a
g=9.8;
Nm = 1E5;
v0 = zeros(Nm,1);
a = zeros(Nm,1);
v0(1) = 1.5 + 0.3*rand(); % v0 in [1.2 1.8]
a(1) = 0.7 + 0.3*rand(); % a in [0.4 1]
%v0(1)=1.3109;
%a(1)=.6116;
DT_theo = repmat([2*v0(1)*a(1) 2*v0(1)*a(1)^2]/g,15,1);
Sigma = 1e-4; % measurement error from correlation?
Pcurr = exp(-1/2/Sigma^2*sum(sum((DT-DT_theo).^2)));
prop_Var_2=1e-5^2; %
sp_2=2.38^2/2; epsilon_2=1e-6;
V_2 = zeros(2,2,Nm);
V_2(:,:,1) = eye(2,2)*prop_Var_2;
R_2 = chol(V_2(:,:,1));
k0_2 = 100;
DT_theo_all=zeros(15,2,Nm);
for k = 2:Nm
    %v0_prop = 1.5 + 0.3*rand();
    %a_prop = 0.7 + 0.3*rand();
    tmp = [v0(k-1);a(k-1)] + (randn(1,2)*R_2)';
    while tmp(1)<1.2 || tmp(1)>1.8 || tmp(2)<0.4 || tmp(2)>1
        tmp = [v0(k-1);a(k-1)] + (randn(1,2)*R_2)';
    end
    v0_prop = tmp(1); a_prop = tmp(2);
    DT_theo = repmat([2*v0_prop*a_prop 2*v0_prop*a_prop^2]/g,15,1);
    Pprop = exp(-1/2/Sigma^2*sum(sum((DT-DT_theo).^2)));
    alpha = min(1,Pprop/Pcurr);
    if (rand()<alpha)
        a(k)=a_prop;
        v0(k)=v0_prop;
        Pcurr = Pprop;
        DT_theo_all(:,:,k)=DT_theo;
    else
        a(k)=a(k-1);
        v0(k)=v0(k-1);
        DT_theo_all(:,:,k)=DT_theo_all(:,:,k-1);
    end
    if mod(k,k0_2)==1
        V_2(:,:,k)=sp_2*cov([v0(1:k-1) a(1:k-1)])+epsilon_2*eye(2,2);
        display(k);
    else
        V_2(:,:,k)=V_2(:,:,k-1);
    end
    R_2=chol(V_2(:,:,k));
end

% plot posterior
burnIn = 1e4;
figure;
subplot(2,1,1)
plot(v0,'k');
set(gca,'fontsize',18,'xticklabel',[])
title('Initial velocity')
ylabel('v_0 (m/s)')
%ylim([1.2 1.6])
subplot(2,1,2)
plot(a,'r')
set(gca,'fontsize',18)
title('Rebound coefficient')
xlabel('Iteration')
ylabel('a')
disp(mean(v0(burnIn:end)))
disp(std(v0(burnIn:end)))
disp(mean(a(burnIn:end)))
disp(std(a(burnIn:end)))
%
figure;
[H,AX,BigAx,P,PAx]=plotmatrix([v0(burnIn:end) a(burnIn:end) ]);
labeltext={'v_0 (m/s)','a'};
for i=1:2
    for j=1:2
        %xlim(AX(i,j),[mu_MT(j)-width,mu_MT(j)+width])
        %ylim(AX(i,j),[mu_MT(i)-width,mu_MT(i)+width])
        set(AX(i,j),'fontsize',18)
    end
    %xlim(PAx(i),[mu_MT(i)-width,mu_MT(i)+width])
    ylabel(AX(i),labeltext{i},'fontsize',18)
    xlabel(AX(i*2),labeltext{i},'fontsize',18)
end
set(gca,'fontsize',18,'fontname','Arial')

%% trajectory
% observation
t=(0:1e-3:.45)';
height = zeros(length(t),9e3);
for k = burnIn+1:Nm
    t1 = v0(k)/g;
    t2 = v0(k)*a(k)/g*2;
    t3 = v0(k)*a(k)^2/g*2;
    h_tmp = [v0(k)^2/2/g-1/2*g*t(t<t1).^2; ...
        v0(k)*a(k)*(t(t>t1&t<t1+t2)-t1)-1/2*g*(t(t>t1&t<t1+t2)-t1).^2; ...
        v0(k)*a(k)^2*(t(t>t1+t2&t<t1+t2+t3)-t1-t2)-1/2*g*(t(t>t1+t2&t<t1+t2+t3)-t1-t2).^2;...
        v0(k)*a(k)^3*(t(t>t1+t2+t3)-t1-t2-t3)-1/2*g*(t(t>t1+t2+t3)-t1-t2-t3).^2];
    height(:,k-burnIn)=h_tmp;
end
%%
t1 = mean(v0(burnIn+1:end))/g;
t2 = mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))/g*2;
t3 = mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2/g*2;
mean_height = [mean(v0(burnIn+1:end))^2/2/g-1/2*g*t(t<t1).^2; ...
    mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))*(t(t>t1&t<t1+t2)-t1)-1/2*g*(t(t>t1&t<t1+t2)-t1).^2; ...
    mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2*(t(t>t1+t2&t<t1+t2+t3)-t1-t2)-1/2*g*(t(t>t1+t2&t<t1+t2+t3)-t1-t2).^2;...
    mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^3*(t(t>t1+t2+t3)-t1-t2-t3)-1/2*g*(t(t>t1+t2+t3)-t1-t2-t3).^2];
%%
figure;
plot(t,height(:,1:100:end),'color',[.5 .5 .5]); hold on
plot(t,mean_height,'r','linewidth',2)
set(gca,'fontsize',18)
xlabel('Time (sec)')
ylabel('Height (m)')
title('Posterior Predicted Trajectory')
ylim([0 .12])
xlim([t(1) t(end)])

% %%
% densBall=8050;
% EBall=180e9;
% NooBall=.305;
% 
% 
% 
% densSamp=4506;
% ESamp=113.8e9;
% NooSamp=0.32;
% 
% 
% RadBall=6.36e-3/2;
% g=9.81;
% 
% 
% t=0:1/fs:1e4/fs;
% 
% 
% v0_bc=[mean(v0(burnIn+1:end)) mean(v0(burnIn+1:end))*mean(a(burnIn+1:end)) mean(v0(burnIn+1:end))*mean(a(burnIn+1:end))^2];
% 
% Force=zeros(size(t,2),3);
% %ForceT=Force';
% DelBall=(1-NooBall^2)/(pi*EBall);
% DelSamp=(1-NooSamp^2)/(pi*ESamp);
% 
% for bc=1:3
%     tc=4.53.*((4.*densBall.*pi.*((DelBall+DelSamp)./3)).^(2./5)).*RadBall.*(v0_bc(bc).^(-1./5));
%     fmax=1.917.*(densBall.^(3./5)).*((DelBall+DelSamp).^(-2./5)).*(RadBall.^2).*v0_bc(bc).^(6./5);  
%     Force(t <= tc,bc)=fmax.*(sin(pi.*t(t <= tc)./tc)).^1.5;
% end
% 
% %%
% figure;
% clr = {'b','r','k'};
% for i = 1:3
%     plot(t*1e6,Force(:,i),clr{i},'linewidth',2); hold on
% end
% set(gca,'fontsize',18)
% xlim([0 30])
% xlabel('Time (\muS)')
% ylabel('Force (N)')
% legend('1st bc','2nd bc','3rd bc')
% 
% %%
% [spec_source,ffc]=plotfft(Force,1/fs);
% for i = 1:3
%     semilogx(ffc/1e3,20*log10(abs(spec_source(:,i))/max(max(abs(spec_source)))),'color',clr{i},'linewidth',2); hold on
% end
% xlim([ffc(1) ffc(end)]/1e3)
% set(gca,'fontsize',16)
% xlabel('Frequency (kHz)')
% ylabel('Spectra Amplitude (dB)')
% legend('1st bc','2nd bc','3rd bc')
% %%
% spec_bc=zeros(8193,16,3);
% TF_bd=zeros(8193,16,3);
% waterlevel = 0;
% for i = 1:3
%     [spec_bc(:,:,i),f]=plotfft(permute(wfm_align_bc1f(:,i,:),[1 3 2]),1/fs);
%     for ch_idx = 1:16
%         TF_bd(:,ch_idx,i)=sqrt((spec_bc(:,ch_idx,i).*conj(spec_source(:,i)))./((abs(spec_source(:,i))).^2 + waterlevel.*repmat(mean(abs(spec_source(:,i)).^2),size(spec_source(:,i),1),1)));
%     end
% end
% %% plot transfer function
% for ch_idx = 1:16
%     if ch_idx<=8
%         %hold on
%         subplot('position',[0.12+l_d(ch_idx)/200 -.05+Axial_d(ch_idx)/80 0.8/10 1.0/5])
%         for i=1:3
%             loglog(f/1e3,abs(TF_bd(:,ch_idx,i)),clr{i},'linewidth',1); hold on
%         end
%         set(gca,'fontsize',14)
%         title(['BD 6, Ch ' num2str(ch_all(ch_idx)) '.'])
%         ylim([.05 1000])
%         xlim([5 1000])
%         grid on
%         %axis off
%     else
%         %hold on
%         subplot('position',[0.12+l(ch_idx-8)/200 -.05+Axial(ch_idx-8)/80 0.8/10 1.0/5])
%         for i=1:3
%             loglog(f/1e3,abs(TF_bd(:,ch_idx,i)),clr{i},'linewidth',1); hold on
%         end
%         set(gca,'fontsize',14)
%         title(['BD 8, Ch ' num2str(ch_all(ch_idx))])
%         ylim([.05 1000])
%         xlim([5 1000])
%         grid on
%         %axis off
%     end
%     if ch_idx == 16 || ch_idx == 6
%         ylabel('Count/N/Hz','fontsize',18)
%     elseif ch_idx == 5 || ch_idx == 13 || ch_idx == 1 || ch_idx == 9
%         xlabel('kHz','fontsize',18)
%     end
% end
% set(gcf,'position',[0 0 1440 960])
% 
% %%
% dtm = 25e-9;
% fnm=1/(2*dtm);
% f1=40e3; f2=500e3;
% [am,bm]=butter(2,[f1/fnm  f2/fnm]);
% figure; plot(0:1/12.5:1/12.5*1e4,wfm_align_bc1f(:,1,8))
% hold on
% plot(0+380:25e-9*1e6:25e-9*(1e4-1)*1e6+380,filtfilt(am,bm,diff(vx1)/25e-9)+1e3)
% plot(0+380:25e-9*1e6:25e-9*(1e4)*1e6+380,filtfilt(am,bm,vx1)*1e6+5e2)
% plot(0+380:25e-9*1e6:25e-9*(1e4)*1e6+380,filtfilt(am,bm,dx1)*1e11-5e2)