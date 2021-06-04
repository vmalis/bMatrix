%GE xml plot

% read all captures
list=dir("*.xml*");
clc
list(1:16,:)=[];
List=list(1:2:end);




data = xml2struct("CaptureWaveform_020.xml.7");




%SSP
%  data.PulseSequence.sequencer{1,1}.Attributes.title

%X
%  data.PulseSequence.sequencer{1,2}.Attributes.title
q=data.PulseSequence.sequencer{1,2}.data.Text;
X=split(q);

X = cellfun(@str2double,X);
X(isnan(X))=[];

timeX = X(1:2:end);
ampX = X(2:2:end);

%Y
%  data.PulseSequence.sequencer{1,3}.Attributes.title
q=data.PulseSequence.sequencer{1,3}.data.Text;
Y=split(q);

Y = cellfun(@str2double,Y);
Y(isnan(Y))=[];

timeY = Y(1:2:end);
ampY = Y(2:2:end);


%Z
%  data.PulseSequence.sequencer{1,4}.Attributes.title
q=data.PulseSequence.sequencer{1,4}.data.Text;
Z=split(q);

Z = cellfun(@str2double,Z);
Z(isnan(Z))=[];

timeZ = Z(1:2:end);
ampZ = Z(2:2:end);


%RHO1
%  data.PulseSequence.sequencer{1,5}.Attributes.title
q=data.PulseSequence.sequencer{1,5}.data.Text;
R=split(q);

R = cellfun(@str2double,R);
R(isnan(R))=[];

timeR = R(1:2:end);
ampR = R(2:2:end);

%RHO2
%  data.PulseSequence.sequencer{1,6}.Attributes.title

%ThETA1
%  data.PulseSequence.sequencer{1,7}.Attributes.title

%ThETA2
%  data.PulseSequence.sequencer{1,8}.Attributes.title



%%%%%% x values
[ampMax,iTimeMax]=findpeaks(ampX);
timeXmax=timeX(iTimeMax);
[ampMin,iTimeMin]=findpeaks(-ampX);
timeXmin=timeX(iTimeMin);
timeXpeaks=cat(1,timeXmax,timeXmin);
ampXpeaks=cat(1,ampMax,-ampMin);

[timeXpeaks,index] = sort(timeXpeaks);
ampXpeaks=ampXpeaks(index);






%%%%%% z values
[ampMax,iTimeMax]=findpeaks(ampZ);
timeZmax=timeZ(iTimeMax);
[ampMin,iTimeMin]=findpeaks(-ampZ);
timeZmin=timeZ(iTimeMin);
timeZpeaks=cat(1,timeZmax,timeZmin);
ampZpeaks=cat(1,ampMax,-ampMin);

[timeZpeaks,index] = sort(timeZpeaks);
ampZpeaks=ampZpeaks(index);




plot(timeZ,ampZ,'r-')
hold on
plot(timeZpeaks,ampZpeaks,'b*')
hold off

figure
plot(timeR,ampR,'r-')
hold off
