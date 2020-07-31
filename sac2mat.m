%% sac to mat
tpt=cell(2,1);
%tst=zeros(8,length(win_good_idx));
%tb=zeros(8,length(win_good_idx));
SacDir='/net/quake/archive/acoustic_emission/AE_mat_balldrop_16ch_@pc/SAC/';
%fids=fopen([SacDir 'sacdatalist'],'r');
for G=1:2
    tpt_tmp=zeros(8,length(tbs));
    for ev=1:length(tbs)
        for ch=1:8
            sacname=['Balldrop.' tdms_name(1:end-5) 'PXI1slot' num2str(board(G)) '_EV' num2str(ev) '.ch' num2str(ch-1) '.SAC'];
            fidd=fopen([SacDir sacname],'r');
            HdrFloats=fread(fidd,70,'float32');
            HdrNhdr=fread(fidd,15,'int32');
            HdrIhdr=fread(fidd,20,'int32');
            HdrLhdr=fread(fidd,5,'int32');
            HeaderStrings=str2mat(fread(fidd,[8 24],'char'));
            B = HdrFloats(6);
            data=fread(fidd,HdrNhdr(10),'float32');
            NZYEAR=HdrNhdr(1);
            NZDAY=HdrNhdr(2);
            NZHOUR=HdrNhdr(3);
            NZMIN=HdrNhdr(4);
            NZSEC=HdrNhdr(5);
            NZMSEC=HdrNhdr(6);
            DELTA=HdrFloats(1);
            TOrigin=HdrFloats(8);
            Atime=HdrFloats(9);
            T0=HdrFloats(11);
            STLA=HdrFloats(32);
            STLO=HdrFloats(33);
            EVLA=HdrFloats(36);
            EVLO=HdrFloats(37);
            EVDP=HdrFloats(39);
            KSTNM=HeaderStrings(1:4,1);
            stanum=str2num(KSTNM(3));
            %matdata(:,ch,idx)=data;
            %tb(ch,idx)=B;
            tpt_tmp(ch,ev)=Atime;
            %tst(ch,idx)=T0;
            fclose(fidd);
        end
    end
    tpt{G}=tpt_tmp;
end