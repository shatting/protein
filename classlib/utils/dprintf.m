function dprintf( s, varargin )
%DPRINTF(s,arg1,arg2,...) = disp(sprintf(s,arg1,arg2,..)

disp(sprintf(s,varargin{:}));

