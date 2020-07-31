% nonlinear leasy square soloution for ball bounces
%  constructo lsq function
x0 = [1.3110; 0.6115; 500e3; 250e3];
lb = [1.3105; 0.6110; 100e3; 10e3];
ub = [1.3115; 0.612; 1000e3; 500e3];
% vn = (v-1.3112)/1e-4; an = (a-1.0.6114)/1e-5; vn = (v-1.3112)/1e-4; vn = (v-1.3112)/1e-4
x = lsqnonlin(@myfun_ball_bounce_alignwfm,x0,lb,ub);
F = myfun_ball_bounce(x);
fs=12.5e6;
dt = 1/fs;
fn=1/(2*dt);
f1=5e3; f2=500e3;
[af,bf]=butter(2,[f1/fn  f2/fn]);
[syn1, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1),x(3), x(4));
[syn2, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1)*x(2),x(3), x(4));
[syn3, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1)*x(2)^2,x(3), x(4));
d_align_syn  = [syn1 syn2 syn3];
% extend to 400 micron, 100 before P and 300 after P
Tp_sample = round(sqrt(sum(r_vec.^2))/Vp*fs/1e3); %samples before P arrival
d_syn_tmp = [zeros(LBefore,3); d_align_syn(Tp_sample:Tp_sample+LAfter,:)];
% d_syn1 = d_syn_tmp(:,1);
% DT_theo = [2*x(1)*x(2) 2*x(1)*x(2)^2]/g;
% DDT_win = cumsum(round(fs*(DT_obs(1,:) - DT_theo(1,:))));
% d_syn2 = [zeros(min(-min(0,DDT_win(1)),LWin),1); d_syn_tmp(max(DDT_win(1)+1,1):min(LWin+DDT_win(1),max(LWin, LWin-DDT_win(1))),2); zeros(min(max(DDT_win(1),0),LWin),1)];
% d_syn3 = [zeros(min(-min(0,DDT_win(2)),LWin),1); d_syn_tmp(max(DDT_win(2)+1,1):min(LWin+DDT_win(2),max(LWin, LWin-DDT_win(2))),3); zeros(min(max(DDT_win(2),0),LWin),1)];
% d_syn = [filtfilt(af,bf,d_syn1); filtfilt(af,bf,d_syn2); filtfilt(af,bf,d_syn3)];
d_syn = reshape(filtfilt(af,bf, d_syn_tmp),LWin*3,1);


%A_obs = max(abs(d_obs)); d_obs_norm = d_obs/A_obs;
%A_syn = max(abs(d_syn)); d_syn_norm = d_syn/A_syn;

d_obs_norm = d_obs/norm(d_obs);
d_syn_norm = d_syn/norm(d_syn);
%%
plot(d_obs_norm,'linewidth',1)
hold on
plot(d_syn_norm,'linewidth',1)