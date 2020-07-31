%%
Dir=pwd;
name={'three','dot_three','fourdot','dot_fourdot','five','dot_five','seven','dot_seven','eight','dot_eight'};
ev_num=zeros(10,1);
ch_most_ev=zeros(10,1);
tbs=cell(10,1);
for G=1:10
    [ev_num(G),ch_most_ev(G)]=max([length(trigs_all{G,1}) length(trigs_all{G,2}) length(trigs_all{G,3}) length(trigs_all{G,4}) length(trigs_all{G,5}) length(trigs_all{G,6}) length(trigs_all{G,7}) length(trigs_all{G,8})]);
    for i=1:ev_num(G)
        if i==1
            tbs{G}=trigs_all{G,ch_most_ev(G)}(1);
        else
            if trigs_all{G,ch_most_ev(G)}(i)-trigs_all{G,ch_most_ev(G)}(i-1)>1e4
                tbs{G}=[tbs{G};trigs_all{G,ch_most_ev(G)}(i)];
            end
        end
    end
end
%%
ev_wfm=cell(10,1);
for G=1:10
    load([Dir '/AE_' name{G} '.mat']);
    fs=10e6;
    for i=1:length(tbs{G})
        Bs=tbs{G}(i)-125e-6*fs;
        Es=tbs{G}(i)+125e-6*fs;
        ev_wfm{G}=[ev_wfm{G} d(Bs:Es,:)];
    end
end