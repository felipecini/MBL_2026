% OD test DAQ input for cheese board (live plot)
clear
clc

%% ===================== OUTPUT DIR ====================
Name = 'test';
OUT_Folder = ['C:\Users\Cardin Lab\Documents\'];
%% ======================================================

% --- DAQ setup ---
dq = daq("ni");
dq.Rate = 2000;
addinput(dq,"cDAQ1Mod1","ai1","Voltage");
addinput(dq,"cDAQ1Mod1","ai0","Voltage");

% --- Figure live ---
All_Data = [];
WinSec = 1;
tBuf = [];
yBuf = [];
y1Buf = [];
fig = figure('Name','LIVE PLOT DAQ - PRESS SPACE TO STOP','NumberTitle','off');

h = animatedline('Color','g','LineWidth',1.5);
h1 = animatedline('Color','b','LineWidth',1.5);
grid on
xlabel("Temps (s)")
ylabel("Voltage (V)")
title("Live plot -- press space bar to stop")


Dat = [];
Dat2 = [];
TS  = [];

t_live = 0;
dt = 1/dq.Rate;


stopFlag = false;

% --- Callback keyboard ---
set(fig,'KeyPressFcn',@(src,event) keyPressCallback(src,event));

start(dq,"continuous")

pause(0.2)  

try
    while ~stopFlag
        data = read(dq, seconds(0.1));
        y = data{:,1};
        y1 = data{:,2};

        Dat = [Dat; y];
        Dat2 = [Dat2; y1];

        t_block = t_live + (1:numel(y))' * dt;
        t_live = t_block(end);
        TS = [TS; t_block];

       % --- append buffer
       tBuf = [tBuf; t_block(:)];
       yBuf = [yBuf; y(:)];
       y1Buf = [y1Buf; y1(:)];

       %--- keep only lax win
       tNow = tBuf(end);
       keep = tBuf >= (tNow - WinSec);

       tBuf = tBuf(keep);
       yBuf =  yBuf(keep);
       y1Buf = y1Buf(keep);

       % --- redraw only in the window
       clearpoints(h);
       clearpoints(h1);
       addpoints(h, tBuf, yBuf);
       addpoints(h1, tBuf, y1Buf);

      xlim([tNow-WinSec, tNow])

       drawnow limitrate
    end
catch ME
    disp("error :")
    disp(ME.message)
end

stop(dq)
All_Data(:,1) = TS;
All_Data(:,2) = Dat; % ai0 -- Mini2p
All_Data(:,3) = Dat2; % ai1 -- Camera

save([OUT_Folder,Name], "All_Data",'-v7.3')
disp(['Stop and save in :', [OUT_Folder,Name],'.mat'])

% --- Fonction callback space bar---
function keyPressCallback(~, event)
    persistent stopFlagHandle
    if isempty(stopFlagHandle)
        stopFlagHandle = evalin('base','stopFlag');
    end
    if strcmp(event.Key,'space')
        assignin('base','stopFlag',true)
    end
end
