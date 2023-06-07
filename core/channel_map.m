

% Channels
time = data(t0,1)- data(1,1); % read time
begin_time = time(1);
end_time = time(end);
Toplid_therm = data(t0,2:7); % top thermistors 0, 60, 120, 180, 240, 300
Botlid_therm = data(t0,8:13); % bottom thermistors 0, 60, 120, 180, 240, 300
Toplid_therm7 = [data(t0,2:7) data(t0,2)];
Botlid_therm7 = [data(t0,8:13) data(t0,8)];
Int_therm = data(t0,14:18); % internal thermistors
intcen_55 = data(t0,14); % center, 55mm
int30_16 = data(t0,15); % 30 degree, 2/3R, 16mm  
int90_105 = data(t0,16); % 90 degree, 2/3R, 105mm 
int210_33 = data(t0,17); % 210 degree, 2/3R, 33mm
int270_16 = data(t0,18); % 270 degree, 2/3R, 16mm
PMU = data(t0,19:20); % Power measurement  
%HallP = data(:,21:22); % 
ExpTank= data(t0,23:25); % Expansion tank
SW_therm = data(t0,26:51); % sidewall ensemble
SW_34 = data(t0,26:31); % 3/4 height
SW_347 = [data(t0,26:31) data(t0,26)];
SW_12 = data(t0,32:43); % 1/2 height
SW_mid13 = [data(t0,32:43) data(t0,32)];
SW_14 = data(t0,44:49); % 1/4 height
SW_147 = [data(t0,44:49) data(t0,44)];
SW_vert0 = [data(t0,50) data(t0,44) data(t0,32) data(t0,26) data(t0,51)]; % 1/8, 1/4, 1/2, 3/4, 7/8
%vert_array = [data(:,8) SW_vert0 data(:,2)];
vert_array = SW_vert0;
SW_loss = [data(t0,32) data(t0,51)]; % sidewall loss in and out
roomtemp = data(t0,53); % room temp