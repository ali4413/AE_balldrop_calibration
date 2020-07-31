% 1. calculate spectra of pieces of signal waveform for some large signals and some small ones
% 2. calculate ratio of spectra of signal of large to signal of small. This should cancel instrument, path etc. and give you source ratio.
% 3. calculate spectra for pieces of noise before or after signals, of same window length
% 4. plot signal/noise and find the frequency range where the signal is at least 3 x noise.
%   (Larger is better. 2x noise means than noise is equal to signal so a bit low. 3 is sort of a good working value)
% 5. Go back to the ratio you calculated in (2) and extract only the part where BOTH signals are good.
%% Test 1: ball drop
wfm_bc_1_test_1 = wfm_bounce(:,:,1);
wfm_bc_2_test_1 = wfm_bounce(:,:,8);
signal_window=[9700 9700+4999];
NL = 5e3; dt = 1/fs;
NFFT=2^nextpow2(NL);
f=(1/(2*dt))*(linspace(0,1,NFFT/2+1));
spec_bc_1_test_1=fft(wfm_bc_1_test_1(signal_window(1):signal_window(2),:),NFFT)/NL;
fft_bc_1_test_1=abs(spec_bc_1_test_1(1:NFFT/2+1,:));
spec_bc_2_test_1=fft(wfm_bc_2_test_1(signal_window(1):signal_window(2),:),NFFT)/NL;
fft_bc_2_test_1=abs(spec_bc_2_test_1(1:NFFT/2+1,:));
loglog(f,fft_bc_2_test_1(:,7))

%%
noise_window=[1 5e3];
wfm_bc_1_test_1_noise = wfm_bc_1_test_1(noise_window(1):noise_window(2),:);
wfm_bc_2_test_1_noise = wfm_bc_2_test_1(noise_window(1):noise_window(2),:);
%NL = noise_window(2)-noise_window(1)+1; dt = 1/fs;
%NFFT=2^nextpow2(NL);
%f=(1/(2*dt))*(linspace(0,1,NFFT/2+1));
spec_bc_1_test_1_noise=fft(wfm_bc_1_test_1_noise,NFFT)/NL;
fft_bc_1_test_1_noise=abs(spec_bc_1_test_1_noise(1:NFFT/2+1,:));
spec_bc_2_test_1_noise=fft(wfm_bc_2_test_1_noise,NFFT)/NL;
fft_bc_2_test_1_noise=abs(spec_bc_2_test_1_noise(1:NFFT/2+1,:));
%%
loglog(f,fft_bc_2_test_1(:,7)); hold on
loglog(f,fft_bc_2_test_1_noise(:,7))
[fcut y]=ginput(2);
%% spectral ratio
fft_ref = fft_bc_1_test_1(:,1);
spectral_ratio_raw = fft_bc_1_test_1./fft_ref;
fn=1/(2*dt);
%................................................
[a,b]=butter(2,[fcut(1)/fn  fcut(2)/fn]);
wfm_bc_1_test_1f=filtfilt(a,b,wfm_bc_1_test_1);
wfm_bc_2_test_1f=filtfilt(a,b,wfm_bc_2_test_1);
spec_bc_1_test_1f=fft(wfm_bc_1_test_1f(signal_window(1):signal_window(2),:),NFFT)/NL;
fft_bc_1_test_1f=abs(spec_bc_1_test_1f(1:NFFT/2+1,:));
spec_bc_2_test_1f=fft(wfm_bc_2_test_1f(signal_window(1):signal_window(2),:),NFFT)/NL;
fft_bc_2_test_1f=abs(spec_bc_2_test_1f(1:NFFT/2+1,:));

fft_ref_f = fft_bc_1_test_1f(:,1);
spectral_ratio_f=fft_bc_1_test_1f./fft_ref_f;
%%
figure;
for ch_idx = 1:16
    if ch_idx<=8
        %hold on
        subplot('position',[0.12+l_d(ch_idx)/200 -.05+Axial_d(ch_idx)/80 0.8/10 1.0/5])
        loglog(f/1e3,fft_bc_1_test_1f(:,ch_idx),'k','linewidth',1)
        set(gca,'fontsize',14)
        title(['BD 6, Ch ' num2str(ch_all(ch_idx)) '.'])
        ylim([.01 100])
        xlim([1 1e3])
        %axis off
    else
        %hold on
        subplot('position',[0.12+l(ch_idx-8)/200 -.05+Axial(ch_idx-8)/80 0.8/10 1.0/5])
        loglog(f/1e3,fft_bc_1_test_1f(:,ch_idx),'r','linewidth',1)
        set(gca,'fontsize',14)
        title(['BD 8, Ch ' num2str(ch_all(ch_idx))])
        ylim([.01 100])
        xlim([1 1e3])
        %axis off
    end
end
set(gcf,'position',[0 0 1440 960])

%%
figure;
for ch_idx = 1:16
    if ch_idx<=8
        %hold on
        subplot('position',[0.12+l_d(ch_idx)/200 -.05+Axial_d(ch_idx)/80 0.8/10 1.0/5])
        loglog(f/1e3,spectral_ratio_f(:,ch_idx),'k','linewidth',1)
        set(gca,'fontsize',18)
        title(['BD 6, Ch ' num2str(ch_all(ch_idx)) '.'])
        ylim([.01 100])
        xlim([1 1e3])
        %axis off
    else
        %hold on
        subplot('position',[0.12+l(ch_idx-8)/200 -.05+Axial(ch_idx-8)/80 0.8/10 1.0/5])
        loglog(f/1e3,spectral_ratio_f(:,ch_idx),'r','linewidth',1)
        set(gca,'fontsize',18)
        title(['BD 8, Ch ' num2str(ch_all(ch_idx))])
        ylim([.01 100])
        xlim([1 1e3])
        %axis off
    end
end
set(gcf,'position',[0 0 1440 960])

%% Test 2: Repeating events from stick-slip exp
ch = 6;
wfm_ref = wfm2(:,ch+1,1);
wfm_all = permute(wfm2(:,ch+1,1:25),[1 3 2]);
signal_window_r=[9700 9700+4999];
fs_r = 25e6;
NL = 5e3; dt = 1/fs_r;
NFFT=2^nextpow2(NL);
f=(1/(2*dt))*(linspace(0,1,NFFT/2+1));
spec_ref=fft(wfm_ref(signal_window(1):signal_window(2),:),NFFT)/NL;
fft_ref=abs(spec_ref(1:NFFT/2+1,:));
spec_all=fft(wfm_all(signal_window(1):signal_window(2),:),NFFT)/NL;
fft_all=abs(spec_all(1:NFFT/2+1,:));
figure; loglog(f,fft_all)

%%
noise_window=[1 5e3];
wfm_all_noise = wfm_all(noise_window(1):noise_window(2),:);
spec_all_noise=fft(wfm_all_noise,NFFT)/NL;
fft_all_noise=abs(spec_all_noise(1:NFFT/2+1,:));

%%
figure;
loglog(f,fft_all./fft_all_noise); 
[fcut_r y]=ginput(2);

%% spectral ratio
fft_ref = fft_all(:,1);
spectral_ratio_raw_r = fft_all./fft_ref;
dt_r = 1/fs_r;
fn_r=1/(2*dt_r);
%................................................
[a,b]=butter(2,[fcut_r(1)/fn_r  fcut_r(2)/fn_r]);
wfm_all_f=filtfilt(a,b,wfm_all);
spec_allf=fft(wfm_all_f(signal_window(1):signal_window(2),:),NFFT)/NL;
fft_allf=abs(spec_allf(1:NFFT/2+1,:));

fft_ref_f_r = fft_allf(:,1);
spectral_ratio_f_r=fft_allf./fft_ref_f_r;
%%
figure; loglog(f,spectral_ratio_raw);
figure; loglog(f,spectral_ratio_f_r);
