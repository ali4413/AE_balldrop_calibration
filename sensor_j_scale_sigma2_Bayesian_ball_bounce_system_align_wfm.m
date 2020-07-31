load Bayes_input.mat
%% ball properties
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
Vs = sqrt(GSamp / densSamp);
%% time range and source-receiver locations
fs=12.5e6;
dt = 1/fs;
t = (0:dt:250e-6)';
rs = [0 0 Length]; % source location
r_vec_all = [xr - rs(1) yr - rs(2) zr - rs(3);
    xr_d - rs(1) yr_d - rs(2) zr_d - rs(3)]; % 16 sensor location
j=2;
r_vec = r_vec_all(j,:);
% MCMC for v0 and a
g=9.8;
Nm = 4E5;
v0 = zeros(Nm,1);
a = zeros(Nm,1);
omega_s = zeros(Nm,1);
varepsilon = zeros(Nm,1);

Ndt = 15*2; % 15 working sensors, 2 bouncing time interval
v0(1) = 1.35;  % v0 in [1.2 1.8] m/s
a(1) = 0.6; % a in [0.4 1]
DT_theo = reshape(repmat([2*v0(1)*a(1) 2*v0(1)*a(1)^2]/g,15,1),Ndt,1);
DT_obs = reshape(DT, Ndt, 1);
omega_s(1) = 364e3;
varepsilon(1) = 23e3;
%omega_s(1) = 550e3 + 450e3*rand(); % omega_s in [100 1000] kHz
%varepsilon(1) = 255e3 + 245e3*rand(); % varepsilon in [10 500] kHz
%v0(1)=1.5;
%a(1)=.55;

d_align=wfm_align_bc1f(:,:,j);
if j == 5
    d_align = -d_align; % sensor 5's polarity is reversed
end
Ppick = Ppick_AIC(j); % P pick with AIC algorithm for 1bc bounce
% cut 400 micron, 100 before P and 300 after P
LWin = round(400e-6*fs)+1;
LBefore = round(100e-6*fs);
LAfter = round(300e-6*fs);
d_obs = reshape(d_align(Ppick - LBefore:round(Ppick + LAfter),:),LWin*3,1);
%d_obs = d_align(round(mean(Pd)) - round(sqrt(sum(r_vec.^2))/Vp*fs/1e3):round(mean(Pd)) + round((400e-6 - sqrt(sum(r_vec.^2))/Vp/1e3)*fs),:);
%
[syn1, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0(1),omega_s(1), varepsilon(1));
[syn2, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0(1)*a(1),omega_s(1), varepsilon(1));
[syn3, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0(1)*a(1)^2,omega_s(1), varepsilon(1));
d_align_syn  = [syn1 syn2 syn3];
% extend to 400 micron, 100 before P and 300 after P
Tp_sample = round(sqrt(sum(r_vec.^2))/Vp*fs/1e3); %samples before P arrival
d_syn_tmp = [zeros(LBefore,3); d_align_syn(Tp_sample:Tp_sample+LAfter,:)];
fn=1/(2*dt);
f1=5e3; f2=500e3;
[af,bf]=butter(2,[f1/fn  f2/fn]);
d_syn = reshape(filtfilt(af,bf, d_syn_tmp),LWin*3,1);
A_obs = max(abs(d_obs)); d_obs_norm = d_obs/norm(d_obs);
A_syn = max(abs(d_syn)); d_syn_norm = d_syn/norm(d_syn);

%%
%nlevDT = 0.0001; % time interval std
%Sigmadt_diag= (DT_obs*nlevDT).^2;
%nlev = 0.2; % 
Sigmadt_diag= (1e-4).^2*ones(Ndt,1);
%Sigma2e_diag=repmat(sum(d_obs_norm.^2)*nlev/(LWin*3),LWin*3,1); %\alpha E/NT
%Sigmadt_diag= [repmat(var(DT_obs(1:15)),15,1);repmat(var(DT_obs(16:30)),15,1)]; % real variance
Sigma2e_diag= repmat(sum(d_obs_norm(1:LBefore).^2)/LBefore,LWin*3,1); % real variance? noise = randn(3*LWin,1)*std(d_obs_norm(1:LBefore));
Sigmadt_inv= diag(1./Sigmadt_diag);
Sigma2e_inv=sparse(diag(1./Sigma2e_diag));
LcurrPrior=-Ndt/2*log(2*pi)-1/2*sum(log(Sigmadt_diag))-1/2*(DT_obs - DT_theo)'*Sigmadt_inv*(DT_obs - DT_theo); % multi-variate Gaussian of Yi^{k}
LcurrLikely = -(LWin*3)/2*log(2*pi)-1/2*sum(log(Sigma2e_diag))-1/2*(d_obs_norm - d_syn_norm)'*Sigma2e_inv*(d_obs_norm - d_syn_norm); % multi-variate Gaussian of Yi^{k}

prop_Var_4=[(2e-3)^2 (1e-3)^2 (1e3)^2 (1e3)^2]; %
sp_4=2.38^2/4; epsilon_4=1e-6;
V_4 = zeros(4,4,Nm);
V_4(:,:,1) = diag(prop_Var_4);
R_4 = chol(V_4(:,:,1));
k0_4 = 10;
tic
%d_sample_all = zeros(3*LWin,Nm);
%d_sample_all(:,1) = d_syn_norm;
for k = 2:Nm
    tmp = [v0(k-1); a(k-1); omega_s(k-1);varepsilon(k-1)] + (randn(1,4)*R_4)';
    while tmp(1)<1.0 || tmp(1)>2 || tmp(2)<0.5 || tmp(2)>0.9||tmp(3)<100e3 || tmp(3)>500e3 || tmp(4)<10e3 || tmp(4)>50e3
        tmp = [v0(k-1); a(k-1); omega_s(k-1);varepsilon(k-1)] + (randn(1,4)*R_4)';
    end
    v0_prop = tmp(1); a_prop = tmp(2);
    omega_s_prop = tmp(3); varepsilon_prop = tmp(4);
    DT_theo_prop = reshape(repmat([2*v0_prop*a_prop 2*v0_prop*a_prop^2]/g,15,1),Ndt,1);
    [syn1_prop, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0_prop,omega_s_prop, varepsilon_prop);
    [syn2_prop, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0_prop*a_prop,omega_s_prop, varepsilon_prop);
    [syn3_prop, ~] = wfm_simulator_fn(densBall,EBall,muBall,densSamp,ESamp,muSamp,RadBall,t,r_vec,v0_prop*a_prop^2,omega_s_prop, varepsilon_prop);
    d_align_syn_prop  = [syn1_prop syn2_prop syn3_prop];
    d_syn_tmp_prop = [zeros(LBefore,3); d_align_syn_prop(Tp_sample:Tp_sample+LAfter,:)];
    d_syn_prop = reshape(filtfilt(af,bf, d_syn_tmp_prop),LWin*3,1);
    %A_syn_prop = max(abs(d_syn_prop));
    d_syn_norm_prop = d_syn_prop/norm(d_syn_prop);
    LPropPrior = -Ndt/2*log(2*pi)-1/2*sum(log(Sigmadt_diag))-1/2*(DT_obs - DT_theo_prop)'*Sigmadt_inv*(DT_obs - DT_theo_prop);
    LPropLikely = -(LWin*3)/2*log(2*pi)-1/2*sum(log(Sigma2e_diag))-1/2*(d_obs_norm - d_syn_norm_prop)'*Sigma2e_inv*(d_obs_norm - d_syn_norm_prop);
    alpha = min(1,exp(LPropLikely + LPropPrior - LcurrLikely - LcurrPrior));
    if (rand()<alpha)
        v0(k)=v0_prop;
        a(k)=a_prop;
        omega_s(k)=omega_s_prop;
        varepsilon(k)=varepsilon_prop;
        LcurrLikely = LPropLikely;
        LcurrPrior = LPropPrior;
        %d_sample_all(:,k) = d_syn_norm_prop;
    else
        v0(k)=v0(k-1);
        a(k)=a(k-1);
        omega_s(k)=omega_s(k-1);
        varepsilon(k)=varepsilon(k-1);
        %d_sample_all(:,k) = d_sample_all(:,k-1);
    end
    if k>1e4
        V_4(:,:,k)=sp_4*cov([v0(1:k-1) a(1:k-1) omega_s(1:k-1) varepsilon(1:k-1)])+epsilon_4*eye(4,4);
        %disp(k);
    else
        V_4(:,:,k)=V_4(:,:,k-1);
    end
    R_4=chol(V_4(:,:,k));
end
%save(['MCMC_4para_scale_sigma_sensor_' num2str(j) '_PpickAIC1bc_test.mat'],'v0','a','omega_s','varepsilon')
toc