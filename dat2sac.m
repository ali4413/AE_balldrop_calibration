function dat2sac(datamat,evname)
% write SAC file from data file: data format wfm x sta
DirSac='/net/quake/archive/acoustic_emission/AE_mat_balldrop_16ch_@pc/SAC/';
chn=8; %channel number
dt=1/12.5e6; %s 
NPTS=length(datamat); % Number of points per data component
%mkdir([DirSac evname]);
for ch=1:chn
    fidsac=fopen([DirSac '/Balldrop.' evname '.ch' num2str(ch-1) '.SAC'],'w');
    HdrFloats=-12345*ones(70,1);
    HdrNhdr=-12345*ones(15,1);
    HdrIhdr=-12345*ones(20,1);
    HdrLhdr=-12345*ones(5,1);
    HeaderString=-12345*ones(24,8);
    NVHDR=6; % Header version number
    B=0; % Beginning time
    E=B+dt*(NPTS-1); % Ending time
    IFTYPE=1; % ITIME{Time series file}
    LEVEN=1; % TRUE, data is evenly spaced
    DELTA=dt; % Time increment
    IZTYPE=11; % Reference time equivalence
    HdrFloats(1)=DELTA;
    HdrFloats(6)=B;
    HdrFloats(7)=E;
    HdrNhdr(7)=NVHDR;
    HdrNhdr(10)=NPTS;
    HdrIhdr(1)=IFTYPE;
    HdrIhdr(3)=IZTYPE;
    HdrLhdr(1)=LEVEN;
    KSTNM=int32(['ch' num2str(ch-1)]);
    KEVNM=int32(evname);
    HeaderString(1:length(KSTNM),1)=KSTNM';
    HeaderString(9:8+length(KEVNM),1)=KEVNM';
    data=datamat(:,ch);
    fwrite(fidsac,HdrFloats,'float32');
    fwrite(fidsac,HdrNhdr,'int32');
    fwrite(fidsac,HdrIhdr,'int32');
    fwrite(fidsac,HdrLhdr,'int32');
    fwrite(fidsac,HeaderString,'char');
    fwrite(fidsac,data,'float32');
    fclose(fidsac);
end