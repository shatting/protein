

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% showtime.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function string=showtime(sec);
% displays time given in seconds (default sec=toc)
% in the format  
%    days:hours:minutes:seconds   
% with 16 characters if <100 days
%
% without output argument, 
% showtime(sec) displays time on screen 
% with output argument, 
% string=showtime(sec) returns time as a string
%
function string=showtime(sec);

if nargin<1, sec=toc; end;

ok=(sec >=0 );

if ~ok,
  % fault in Matlab's etime routine
  % month boundary passed
  str='new month; reset';
  tic;
end;

  

while ok,
  if sec>=60,
    min=fix(sec/60);
    sec=fix(sec-60*min+0.5);
    if min>=60,
      h=fix(min/60);
      min=min-60*h;
      if h>=24,
        d=fix(h/24);
        h=h-24*d;
        str=sprintf('%2i:%2.2i:%2.2i:%2.2i days',d,h,min,sec);
        break;
      end;
      str=sprintf('  %2i:%2.2i:%2.2i hours',h,min,sec);
      break;
    end;
    str=sprintf('%5i:%2.2i minutes',min,sec);
    break;
  end;
  str=sprintf('%12.3f sec',sec);
  break;
end;

if nargout==0,
  disp(str)
else
  string=str;
end;

