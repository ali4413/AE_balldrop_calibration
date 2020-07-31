for G=1:2
    for ev=1:length(tbs)
        datamat=ev_wfm{G,ev};
        evname=[tdms_name(1:end-5) 'PXI1slot' num2str(board(G)) '_EV' num2str(ev)];
        dat2sac(datamat,evname);
        fclose all;
    end
end