function [amp, time, a,t,iTime] = gradTimePoints(data, name)
%==========================================================================
% subroutine to mark all the time points with gradient changes
% 
%   input:  structure from GE xml file
%           gradient axis
%
%   output: all the amplitude and time from xml
%           important gradient switch points
%
%==========================================================================
% 03/2019 - VM (vmalis@ucsd.edu)   
%==========================================================================


switch name
    case 'X'
        idx = 2;
    case 'Y'
        idx = 3;
    case 'Z'
        idx = 4;
    otherwise
        return
end

q=data.PulseSequence.sequencer{1,idx}.data.Text;
G=split(q);

G = cellfun(@str2double,G);
G(isnan(G))=[];

time = G(1:2:end);
time(end)=[];
ampG = G(2:2:end);
ampG(end)=[];
amp = ampG;

[ampMax,iTimeMax]=findpeaks(ampG);
ampMax=cat(1,ampMax,ampG(iTimeMax+1));
timePeakMax=time(cat(1,iTimeMax,iTimeMax+1));

[ampMin,iTimeMin]=findpeaks(-ampG);
ampMin=cat(1,-ampMin,ampG(iTimeMin+1));
timePeakMin=time(cat(1,iTimeMin,iTimeMin+1));

[iTimeZ,~]=find(~ampG);
timeZ = time(iTimeZ);

iTime = cat(1,iTimeMax,iTimeMax+1,iTimeMin,iTimeMin+1,iTimeZ);

timePeaks=cat(1,timePeakMax,timePeakMin,timeZ);
ampPeaks=cat(1,ampMax,ampMin,ampG(iTimeZ));

[timePeaks,index] = sort(timePeaks);
ampPeaks=ampPeaks(index);
iTime=iTime(index);
[iTime,unq]=unique(iTime);

a=ampPeaks(unq);
t=timePeaks(unq);

a(1)=[];
t(1)=[];
iTime(1)=[];

end