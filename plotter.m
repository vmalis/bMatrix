%GE xml plot

% read all captures
list=dir("*.xml*");
nGrFiles = size(list,1);

f = waitbar(0,'Reading gradient data...');

for i=1:size(list,1)
    waitbar(i/nGrFiles,f,'Reading gradient data...');
    GradData(i).data     = xml2struct(list(i).name);
    GradData(i).Delta    = str2num(list(i).name(17:19));
    GradData(i).GrNumber = str2num(list(i).name(25))-2;
   
    data = GradData(i).data;
    % X - read out
    [amp, time, a,t,iTime] = gradTimePoints(data,'X');
    GradData(i).xAmp   = amp;
    GradData(i).xTime  = time;
    GradData(i).xAp    = a;
    GradData(i).xt     = t;
    GradData(i).iTimeX = iTime;
    
    % Y - phase encode
    [amp, time, a,t,iTime] = gradTimePoints(data,'Y');
    GradData(i).yAmp   = amp;
    GradData(i).yTime  = time;
    GradData(i).yAp    = a;
    GradData(i).yt     = t;
    GradData(i).iTimeY = iTime;
    
    % Z - phase encode
    [amp, time, a,t,iTime] = gradTimePoints(data,'Z');
    GradData(i).zAmp   = amp;
    GradData(i).zTime  = time;
    GradData(i).zAp    = a;
    GradData(i).zt     = t;
    GradData(i).iTimeZ = iTime;
    
end

close(f)


%% plots

N=1;

figure

subplot(3,1,1)
plot(GradData(N).xTime,GradData(N).xAmp,'r-')
hold on
plot(GradData(N).xt,GradData(N).xAp,'b*')
hold off

subplot(3,1,2)
plot(GradData(N).yTime,GradData(N).yAmp,'r-')
hold on
plot(GradData(N).yt,GradData(N).yAp,'b*')
hold off

subplot(3,1,3)
plot(GradData(N).zTime,GradData(N).zAmp,'r-')
hold on
plot(GradData(N).zt,GradData(N).zAp,'b*')
hold off



%% get parameters and save to structure

for n=1:size(GradData,2)


    % READOUT--------------------------------------

    if size(GradData(n).xAp,1)>160  %with diffusion

        Gradients(n).Gdr  = GradData(n).xAp(2);       
        Gradients(n).Gcr  = GradData(n).xAp(6);
        Gradients(n).Grdp = GradData(n).xAp(22);
        Gradients(n).Gro  = GradData(n).xAp(30);

    else                            %no diffusion

        Gradients(n).Gdr  = 0;       
        Gradients(n).Gcr  = GradData(n).xAp(2);
        Gradients(n).Grdp = GradData(n).xAp(14);
        Gradients(n).Gro  = GradData(n).xAp(22);

    end

    % PhaseEncode------------------------------------

    if size(GradData(n).yAp,1)>160  %with diffusion

        Gradients(n).Gdp  = GradData(n).yAp(2);       
        Gradients(n).Gcp  = GradData(n).yAp(6);
        Gradients(n).Gpdp = GradData(n).yAp(22);
        Gradients(n).Gpe  = GradData(n).yAp(26);

    else                            %no diffusion

        Gradients(n).Gdp  = 0;       
        Gradients(n).Gcp  = GradData(n).yAp(2);
        Gradients(n).Gpdp = GradData(n).yAp(14);
        Gradients(n).Gpe  = GradData(n).yAp(18);

    end

    % SliceSelect------------------------------------

    if size(GradData(n).zAp,1)>88  %with diffusion

        Gradients(n).Gsl  = GradData(n).zAp(6);
        Gradients(n).Grf  = GradData(n).zAp(59);
        Gradients(n).Gds  = GradData(n).zAp(62);
        Gradients(n).Gcs  = GradData(n).zAp(66);
        Gradients(n).Gsl2 = GradData(n).zAp(82);

    else                            %no diffusion

        Gradients(n).Gsl  = GradData(n).zAp(6);
        Gradients(n).Grf  = GradData(n).zAp(59);
        Gradients(n).Gds  = 0;
        Gradients(n).Gcs  = GradData(n).zAp(62);
        Gradients(n).Gsl2 = GradData(n).zAp(78);

    end
    
    
    % t
    
    % same for all
    Gradients(n).TE = (GradData(n).xt(end-71)-(GradData(n).xt(end-71) - ...
    GradData(n).xt(end-72))/2) - (GradData(n).zt(60)-GradData(n).zt(1))/2;
    Gradients(n).t5s  = GradData(n).zt(57); 
    Gradients(n).d1   =  GradData(n).zt(3)-GradData(n).zt(1);
    Gradients(n).d5s  =  GradData(n).zt(59)-GradData(n).zt(57);
    Gradients(n).d7   =  GradData(n).xt(end-1)-GradData(n).xt(end-3);
    Gradients(n).eps4 =  0;       
    Gradients(n).eps5s  =  GradData(n).zt(58)-GradData(n).zt(57);
    Gradients(n).eps6   =  GradData(n).yt(end-7)-GradData(n).yt(end-8);
    Gradients(n).eps7   =  GradData(n).xt(end-2)-GradData(n).xt(end-3);
    
    Gradients(n).t6 = GradData(n).yt(end-7:-4:end-4*36);
    
    if size(GradData(n).xAp,1)>160  %with diffusion
        Gradients(n).t21  = GradData(n).xt(1);
        Gradients(n).t22  = GradData(n).xt(17); 
        Gradients(n).eps2   =  GradData(n).xt(3)-GradData(n).xt(1);
        Gradients(n).t31  = GradData(n).xt(5);
        Gradients(n).t32  = GradData(n).xt(13);
        Gradients(n).t41  = GradData(n).xt(8);
        Gradients(n).t5rp = GradData(n).xt(21);
        Gradients(n).t71  = GradData(n).xt(25);
        Gradients(n).d2   =  GradData(n).xt(3)-GradData(n).xt(1);
        Gradients(n).d3   =  GradData(n).xt(7)-GradData(n).xt(5);
        Gradients(n).d4   =  GradData(n).xt(9)-GradData(n).xt(8);
        Gradients(n).d5rp =  GradData(n).xt(23)-GradData(n).xt(21);
        Gradients(n).eps3   =  GradData(n).xt(6)-GradData(n).xt(5);
        Gradients(n).eps5rp =  GradData(n).xt(22)-GradData(n).xt(21);
    else
        Gradients(n).t21  = 0;
        Gradients(n).t22  = 0;
        Gradients(n).eps2 = 0;
        Gradients(n).t31  = GradData(n).xt(1);
        Gradients(n).t32  = GradData(n).xt(9);
        Gradients(n).t41  = GradData(n).xt(4);
        Gradients(n).t5rp = GradData(n).xt(13);
        Gradients(n).t71  = GradData(n).xt(17); 
        Gradients(n).d2 = 0;
        Gradients(n).d3 = GradData(n).xt(3)-GradData(n).xt(1);
        Gradients(n).d4 = GradData(n).xt(5)-GradData(n).xt(4);
        Gradients(n).d5rp = GradData(n).xt(15)-GradData(n).xt(13);
        Gradients(n).eps3 = GradData(n).xt(2)-GradData(n).xt(1);
        Gradients(n).eps5rp = GradData(n).xt(14)-GradData(n).xt(13);
    end
    
    
    if size(GradData(n).zAp,1)>88 %with diffusion
        Gradients(n).t42  = GradData(n).zt(end-14);   
    else
        Gradients(n).t42  = GradData(n).zt(end-10);
    end
    
end

b = zeros(3,3,size(GradData,2));

for n=1:size(GradData,2)
     b(:,:,n)= b_matrix(Gradients(n));
end