function [ s ] = ask( q, t )
%ASK Summary of this function goes here
%   Detailed explanation goes here

if nargin==2 && t > 0,

    if 0,
        try 
            dprintf('%s [%i seconds timeout, ctrl-c for NO, else YES',q,t);
            pause(t);
            s = true;
            return
        catch
            s = false;
            return
        end
    else
        dprintf('%s [waiting %i secs. {0} = NO, {any other or timeout} = YES]',q,t);
        k = getkeywait(t);
        if k > 0 && char(k)=='0',
            s = false;
        else
            s = true;
        end
    end
else
    s = input([q,' [0=NO, any other=YES]: '],'s');
    if (strcmp(s,'0')),
        s = false;
    else
        s = true;
    end
end