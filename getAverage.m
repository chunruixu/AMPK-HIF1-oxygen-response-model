function [Ave, steadystate, oscillator] = getAverage(yout_, tout,timeBegin, timeEnd)
% yout_: specific species; tout: time
% whether oscillation or steady state
steadystate = 0; oscillator = 0;
%  if isempty(timeBegin)
%      timeBegin = 50;
%  end
%  if isempty(timeEnd)
%      timeEnd = tout(end);
%  end
 
 A = find(tout>=timeBegin);
 timeBegin = tout(A(1));
 A = find(tout==timeBegin);%index of begining time
 
 B = find(tout<=timeEnd);
 timeEnd = tout(B(end));
 B = find(tout==timeEnd);%index of ending time
 y = yout_(A:B);
 t = tout(A:B);
 [pks,locs_pks] = findpeaks(y);%check peaks
 lmin = islocalmin(y);
 locs_lmin = find(lmin==1);
 lmin = y(locs_lmin);
 All = [pks' lmin'];
 
 if std(y)<=0.1 || std(All)<=0.1 %if the std is less/equal to 10%, the sim is regarded as achieving the steady state
     steadystate = 1;
     Ave = mean(y);%average
 elseif std(pks)<=0.1 && std(lmin)<=0.1
     oscillator = 1;
     Ave = mean(y(locs_pks(1):locs_pks(end)));
 else
     Ave = mean(y);
     warning('sim does not arrive steady-state or oscillation')
 end
 
%  if steadystate == 1
%      Ave = mean(y);
%  elseif oscillator == 1
%      Ave = mean(y(locs_pks(1):locs_pks(end)));
%  else
% %      Ave = 'NaN';
%      Ave = mean(y);
%      warning('sim does not arrive steady-state or oscillation')
%  end
end
 