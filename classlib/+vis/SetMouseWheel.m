function setMouseWheel(fcn,HWIN) 
% setMouseWheel - set callback for mouseWheel for figure 
% 
% INPUT 
% fcn: calbback function 
% HWIN: handle of figure 
% 
%TESTED WITH 
% matlab 7.0.1 
% 
% Nanne van der Zijpp 
% Modelit 
% www.modelit.nl 


if nargin<1 
    fcn=@testcallback; 
end 
if nargin<2 
    HWIN=gcf; 
end 


%get RootPane for specified fure 
RootPane=getRootPane(HWIN); 
setappdata(HWIN,'RootPane',RootPane); %bind RootPane to prevent if from being destroyed 
set(RootPane,'MouseWheelMovedCallback',fcn); 
%-------------------- 
function testcallback(varargin) 
disp(sprintf('testcallback called with %d input arguments',length(varargin))); 
for k=1:length(varargin) 
    disp(sprintf('argument %d:',k)); 
    disp(varargin{k}); 
end 
%-------------------- 
function RootPane=getRootPane(HWIN) 
%get RootPane for specified figure 


if nargin<1 
    HWIN=gcf; 
end 
%create invisible dummy java object and wrap it in a java component 
jobj=javax.swing.JLabel; 
jobj.setOpaque(0); 


[jobj,h] = javacomponent(jobj,[],HWIN); 
RootPane=[]; 
T0=now; 
while isempty(RootPane) 
    pause(.01); 
    RootPane=jobj.getRootPane; 
    if (now-T0)*60*1440 > 1 
        error('Time out while retreiving root pane'); 
    end 
end 
delete(h); 


