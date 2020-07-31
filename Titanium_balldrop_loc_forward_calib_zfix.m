%% Solnhofen sensor info
Rho = 4500; %g/cm^3
Vp = 6560; % m/s
%Vs = 3120; % m/s

%% sensor info
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
xr = Diameter/2 * sin(Azimuth);
yr = Diameter/2 * cos(Azimuth);
zr = Axial;
xr_d = Diameter/2 * sin(Azimuth_d);
yr_d = Diameter/2 * cos(Azimuth_d);
zr_d = Axial_d;



%%
idx=0;
%loc_grids = zeros(3,83481021);
for xs = -ceil(Diameter/2):1:ceil(Diameter/2)
    for ys = -ceil(Diameter/2):1:ceil(Diameter/2)
        for zs = Length
            %if xs^2+ys^2<=(Diameter/2)^2
                idx = idx+1;
                %loc_grids(:,idx) = [xs; ys; zs];
            %end
        end
    end
end
num_grids=idx;
%%
travel_time_grids_P = zeros(8,num_grids);
%travel_time_grids_S = zeros(8,num_grids);
loc_grids = zeros(3,idx);
idx = 0;
%%
for xs = -ceil(Diameter/2):1:ceil(Diameter/2)
    for ys = -ceil(Diameter/2):1:ceil(Diameter/2)
        for zs = Length
            %if xs^2+ys^2<=(Diameter/2)^2
                idx = idx +1; 
                loc_grids(:,idx) = [xs; ys; zs];
                travel_time_grids_P(:,idx) = sqrt((xr-xs).^2 + (yr-ys).^2 + (zr-zs).^2)/Vp; % ms
                %travel_time_grids_S(:,idx) = sqrt((xr-xs).^2 + (yr-ys).^2 + (zr-zs).^2)/Vs; % ms
            %end
        end
    end
end

%%
travel_time_grids_P_d = zeros(8,num_grids);
%travel_time_grids_S_d = zeros(8,num_grids);
loc_grids = zeros(3,num_grids);
idx = 0;
%%
for xs = -ceil(Diameter/2):1:ceil(Diameter/2)
    for ys = -ceil(Diameter/2):1:ceil(Diameter/2)
        for zs = Length
            %if xs^2+ys^2<=(Diameter/2)^2
                idx = idx +1; 
                loc_grids(:,idx) = [xs; ys; zs];
                travel_time_grids_P_d(:,idx) = sqrt((xr_d-xs).^2 + (yr_d-ys).^2 + (zr_d-zs).^2)/Vp; % ms
                %travel_time_grids_S_d(:,idx) = sqrt((xr_d-xs).^2 + (yr_d-ys).^2 + (zr_d-zs).^2)/Vs; % ms
            %end
        end
    end
end

