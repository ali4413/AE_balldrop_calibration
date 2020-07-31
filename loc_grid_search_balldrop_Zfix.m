%tpt_theo=zeros(8,96);
good_loc_idx=cell(10,1);
t0_all=cell(10,1);
%good_loc_idx_P=zeros(42,1);
for G=1:10
    loc_idx_tmp=zeros(length(tbs{G}),1);
    t0_tmp=zeros(length(tbs{G}),1);
    for ev=1:length(tbs{G})
        t0=zeros(num_grids,1);
        error2=zeros(num_grids,1);
        error2_t0=zeros(num_grids,1);
        %vp=1500; %m/s
        minerror=1;
        %tst_temp=tst(:,ev);
        tpt_temp=tpt{G}(:,ev);
        for loc_idx=1:num_grids
            if mod(G,2)==1
                %t0(loc_idx)=mean(tst_temp(tst_temp~=-12345)-travel_time_grids_S(tst_temp~=-12345,loc_idx)*1e-3);
                %error2_t0(loc_idx)=norm(travel_time_grids_S(tst_temp~=-12345,loc_idx)*1e-3+t0(loc_idx)-tst_temp(tst_temp~=-12345));
                %error2(loc_idx)=norm(travel_time_grids_S(tst_temp~=-12345,loc_idx)*1e-3-tst_temp(tst_temp~=-12345));
                t0(loc_idx)=mean(tpt_temp(tpt_temp~=-12345&tpt_temp~=0)-travel_time_grids_P(tpt_temp~=-12345&tpt_temp~=0,loc_idx)*1e-3);
                error2_t0(loc_idx)=norm(travel_time_grids_P(tpt_temp~=-12345&tpt_temp~=0,loc_idx)*1e-3+t0(loc_idx)-tpt_temp(tpt_temp~=-12345&tpt_temp~=0));
                error2(loc_idx)=norm(travel_time_grids_P(tpt_temp~=-12345&tpt_temp~=0,loc_idx)*1e-3-tpt_temp(tpt_temp~=-12345&tpt_temp~=0));
                if error2_t0(loc_idx)<minerror
                    minerror=error2_t0(loc_idx);
                    loc_opt_grid=loc_idx;
                end
            else
                t0(loc_idx)=mean(tpt_temp(tpt_temp~=-12345&tpt_temp~=0)-travel_time_grids_P_d(tpt_temp~=-12345&tpt_temp~=0,loc_idx)*1e-3);
                error2_t0(loc_idx)=norm(travel_time_grids_P_d(tpt_temp~=-12345&tpt_temp~=0,loc_idx)*1e-3+t0(loc_idx)-tpt_temp(tpt_temp~=-12345&tpt_temp~=0));
                error2(loc_idx)=norm(travel_time_grids_P_d(tpt_temp~=-12345&tpt_temp~=0,loc_idx)*1e-3-tpt_temp(tpt_temp~=-12345&tpt_temp~=0));
                if error2_t0(loc_idx)<minerror
                    minerror=error2_t0(loc_idx);
                    loc_opt_grid=loc_idx;
                end
            end
        end
        %         figure
        %         %plot(travel_time_grids_S(tst_temp~=-12345,loc_opt_grid)*1e-3+t0(loc_opt_grid),'-o'); hold on
        %         %plot(tst_temp(tst_temp~=-12345),'-ro')
        %         plot(travel_time_grids_P(tpt_temp~=-12345&tpt_temp~=0,loc_opt_grid)*1e-3+t0(loc_opt_grid),'-o'); hold on
        %         plot(tpt_temp(tpt_temp~=-12345&tpt_temp~=0),'-ro')
        %
        %        set(gca,'fontsize',16)
        %        xlabel('# of P picks')
        %        ylabel('P arrival time')
        loc_idx_tmp(ev)=loc_opt_grid;
        t0_tmp(ev)=t0(loc_opt_grid);
        %good_loc_idx_P(ev)=loc_opt_grid;
    end
    good_loc_idx{G}=loc_idx_tmp;
    t0_all{G}=t0_tmp;
end
% %%
% figure
% for i=1:8
%     subplot(4,2,i)
%     tpt_temp=tpt(:,i);
%     loc_opt_grid=good_loc_idx(i);
%     plot(travel_time_grids_P(tpt_temp~=-12345&tpt_temp~=0,loc_opt_grid)*1e-3+t0_all(i),'-o'); hold on
%     plot(tpt_temp(tpt_temp~=-12345&tpt_temp~=0),'-ro')
%     set(gca,'fontsize',16)
%     if i==7 || i==8
%         xlabel('# of P picks')
%     else
%         set(gca,'xticklabel',[])
%     end
%     title(['source: ch' num2str(i-1)])
%     ylabel('P arrival time')
% end
%% plot cylinder
in2mm=25.4;
Diameter = 1.5*in2mm+8; % mm
Length = 3*in2mm; % mm
%sensor_id = 1:8;
Axial = [12/5; 9/5; 6/5; 3/5; 12/5; 9/5; 6/5; 3/5]*in2mm; % mm
Azimuth = [pi/2; pi/4; 0; 7*pi/4; 3*pi/2; 5*pi/4; pi; 3*pi/4];
Axial_d = [12/5; 9/5; 6/5; 3/5; 12/5; 9/5; 6/5; 3/5]*in2mm; % mm
Azimuth_d = [0; 7*pi/4; 3*pi/2; 5*pi/4; pi; 3*pi/4; pi/2; pi/4];
sensor_id = 8:-1:1;
channel_id = 0:7;

figure
R=[Diameter/2 Diameter/2];
x0=0; y0=0; z0=0;
N=100;
[X,Y,Z] = cylinder(R,N);
testsubjectO = surf(X+x0,Y+y0,Z*Length+z0);
set(testsubjectO,'FaceAlpha',0.2,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)
axis equal
set(gca,'fontname','arial','fontsize',16)
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
hold on
% plot sensor
% plot3(xr, yr, zr, 'ro', 'markersize', 10, 'markerfacecolor','r','markeredgecolor','r')
% for i = 1:8
%     text(xr(i)+ 2.5, yr(i), zr(i),[num2str(channel_id(i)) '(' num2str(sensor_id(i)) ')'],'fontsize',16)
% end
%
% plot sources
%mis_loc=zeros(8,3);
color={'r','g','b','k','m'};
for G=1:10
    %for ev=1:length(tbs{G})
    for ev=1
        if mod(G,2)==1
            plot3(loc_grids(1,good_loc_idx{G}(ev)),loc_grids(2,good_loc_idx{G}(ev)), loc_grids(3,good_loc_idx{G}(ev)), 'p', 'markersize', 10, 'markerfacecolor',color{(G+1)/2},'markeredgecolor',color{(G+1)/2})
        else
            plot3(loc_grids(1,good_loc_idx{G}(ev)),loc_grids(2,good_loc_idx{G}(ev)), loc_grids(3,good_loc_idx{G}(ev)), '^', 'markersize', 10, 'markerfacecolor',color{G/2},'markeredgecolor',color{G/2})
        end
        %mArrow3([xr(ev) yr(ev) zr(ev)],loc_grids(:,good_loc_idx(ev))', 'linewidth',2, 'color', 'black')
        %mis_loc(ev,:)=abs([xr(ev)-loc_grids(1,good_loc_idx(ev)) yr(ev)-loc_grids(2,good_loc_idx(ev)) zr(ev)-loc_grids(3,good_loc_idx(ev))]);
    end
end
hold on
r1=Diameter/2; %inch
theta1=0:pi/18:2*pi;
x1=r1*cos(theta1); y1=r1*sin(theta1);
testsubject3 = patch(x1+x0,y1+y0,12*ones(size(x1)),'b');
set(testsubject3,'FaceAlpha',0.2,'EdgeColor','red','EdgeAlpha',0,'DiffuseStrength',1,'AmbientStrength',1)

