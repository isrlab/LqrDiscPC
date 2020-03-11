function [sys, x0, str,ts] = f16sfcn(t, x ,u, flag,x0)
switch flag
    case 0 % initialize
        str =[];
        ts = [0 0];           
        s = simsizes ;   
        s.NumContStates = 4;
        s.NumDiscStates = 0;
        s.NumOutputs = 4;
        s.NumInputs = 2;
        s.DirFeedthrough = 0;
        s.NumSampleTimes = 1; 
        sys = simsizes(s);
 
    case 1  
        sys = f16longi(x,u);        
         
    case 3  % output
        sys = x ;
        
    otherwise
        sys = [];
end