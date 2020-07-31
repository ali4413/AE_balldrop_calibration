% nonlinear leasy square soloution for ball bounces
%  constructo lsq function
omega_opt = 264e3;
varepsilon_opt = 86e4;
x0 = [1.3112; 0.6113];
x = lsqnonlin(@myfun_ball_bounce_v0_a,x0);
F = myfun_ball_bounce_v0_a(x);
[syn1, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1),omega_opt,varepsilon_opt);
[syn2, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1)*x(2),omega_opt,varepsilon_opt);
[syn3, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1)*x(2)^2,omega_opt,varepsilon_opt);
d_align_syn  = [syn1 syn2 syn3];
% extend to 400 micron, 100 before P and 300 after P
Tp_sample = round(sqrt(sum(r_vec.^2))/Vp*fs/1e3); %samples before P arrival
d_syn_tmp = [zeros(LBefore,3); d_align_syn(Tp_sample:Tp_sample+LAfter,:)];
d_syn1 = d_syn_tmp(:,1);
DT_theo = [2*x(1)*x(2) 2*x(1)*x(2)^2]/g;
DDT_win = cumsum(round(fs*(DT_obs(1,:) - DT_theo(1,:))));
d_syn2 = [zeros(min(-min(0,DDT_win(1)),LWin),1); d_syn_tmp(max(DDT_win(1)+1,1):min(LWin+DDT_win(1),max(LWin, LWin-DDT_win(1))),2); zeros(min(max(DDT_win(1),0),LWin),1)];
d_syn3 = [zeros(min(-min(0,DDT_win(2)),LWin),1); d_syn_tmp(max(DDT_win(2)+1,1):min(LWin+DDT_win(2),max(LWin, LWin-DDT_win(2))),3); zeros(min(max(DDT_win(2),0),LWin),1)];
d_syn = [d_syn1; d_syn2; d_syn3];

A_obs = max(abs(d_obs)); d_obs_norm = d_obs/A_obs;
A_syn = max(abs(d_syn)); d_syn_norm = d_syn/A_syn;

plot(d_obs_norm)
hold on
plot(d_syn_norm)