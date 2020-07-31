function F = myfun_ball_bounce_v0_a(x)
load obs.mat d_obs DT_obs Length xr yr zr xr_d yr_d zr_d rs
omega_opt = 264e3;
varepsilon_opt = 86e4;

densBall = 8050; % Density
EBall = 180e9; % Young's modulas
muBall = .305; % Possion's ratio

densSamp = 4506; % Density
ESamp = 113.8e9; % Young's modulas
muSamp = 0.32; % Possion's ratio

KSamp = ESamp / (3 * (1 - 2 * muSamp));
GSamp = ESamp / (2 * (1 + muSamp));

RadBall = 6.36e-3/2;
%Vp = 6560; % m/s; Uli measured; Vp = 6070 from web as Vs
%Vs = 3310; % m/s; https://www.nde-ed.org/GeneralResources/MaterialProperties/UT/ut_matlprop_metals.htm
Vp = sqrt((KSamp + 4 * GSamp/ 3) / densSamp);

fs=12.5e6;
dt = 1/fs;
t = (0:dt:250e-6)';
rs = [0 0 Length]; % source location
r_vec_all = [xr - rs(1) yr - rs(2) zr - rs(3);
    xr_d - rs(1) yr_d - rs(2) zr_d - rs(3)]; % 16 sensor location
j=4;
r_vec = r_vec_all(j,:);
g=9.8;
% cut 400 micron, 100 before P and 300 after P
LWin = round(400e-6*fs)+1;
LBefore = round(100e-6*fs);
LAfter = round(300e-6*fs);
[syn1, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1),omega_opt, varepsilon_opt);
[syn2, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1)*x(2),omega_opt, varepsilon_opt);
[syn3, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,x(1)*x(2)^2,omega_opt, varepsilon_opt);
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

F = d_obs_norm-d_syn_norm;
