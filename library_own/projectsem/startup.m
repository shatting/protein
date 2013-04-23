if 0,
  dir='C:\Users\browny\Documents\uni\da\users\neum\matlab\';
else
  dir='/users/neum/matlab/';
end;

rand('state',sum(100*clock));
addpath([dir,'ai/suffclass']);
addpath([dir,'util']) 
addpath([dir,'stat']) 
addpath([dir,'protein/basic']) 
%addpath([dir,'intlab']) 