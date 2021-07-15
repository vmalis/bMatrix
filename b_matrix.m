%==========================================================================
%% Calculate full b-matrix for the given gradient shapes
%==========================================================================
%   06/2021 - VM (vmalis@ucsd.edu)   
%==========================================================================
function [b,D2] = b_matrix(Gradients)


%  b-values common units are [s/mm^2]
%  GE gives:
%       time        [us]
%       gradients   [G/cm]

%  Conversion factors:
%       time        -> [us]        = 1E-6 [s]               =>   1E-6
%       gradients   -> [G/cm]      = 1E-4 [T] / 1E1 [mm]    =>   1E-5
%       gamma       -> [2pi*MHz/T] = 2pi*1E6/(T*s)          =>   2*pi*1E6


%--------------------------------------------------------------------------
%% INPUT
res   = 35;
gamma = 2*pi*1E6*42.577478518;
gsq   = gamma^2;



% gradients  [G/cm] = 1E-4 T / 1E-2 m = 1E-2 T/m = 
Gdr  = Gradients.Gdr*1E-5;       
Gcr  = Gradients.Gcr*1E-5;
Grdp = Gradients.Grdp*1E-5;
Gro  = Gradients.Gro*1E-5;
Gdp  = Gradients.Gdp*1E-5;
Gcp  = Gradients.Gcp*1E-5;
Gpdp = Gradients.Gpdp*1E-5;
Gpe  = Gradients.Gpe*1E-5;
Gsl  = Gradients.Gsl*1E-5;
Gds  = Gradients.Gds*1E-5;
Gcs  = Gradients.Gcs*1E-5;
Gsl2 = Gradients.Gsl2*1E-5;
Grf  = Gradients.Grf*1E-5;

% t
TE   = Gradients.TE*1E-6;
t21  = Gradients.t21*1E-6;
t22  = Gradients.t22*1E-6;
t31  = Gradients.t31*1E-6;
t32  = Gradients.t32*1E-6;
t41  = Gradients.t41*1E-6;
t42  = Gradients.t42*1E-6;
t5rp = Gradients.t5rp*1E-6;
t5s  = Gradients.t5s*1E-6;
t71  = Gradients.t71*1E-6;

% delta
d1   =  Gradients.d1*1E-6;
d2   =  Gradients.d2*1E-6;
d3   =  Gradients.d3*1E-6;
d4   =  Gradients.d4*1E-6;
d5rp =  Gradients.d5rp*1E-6;
d5s  =  Gradients.d5s*1E-6;
d7   =  Gradients.d7*1E-6;
%d_sl = 3200*1E-6;% slice select

% eps
eps2   =  Gradients.eps2*1E-6;
eps3   =  Gradients.eps3*1E-6;
eps4   =  Gradients.eps4*1E-6;
eps5rp =  Gradients.eps5rp*1E-6;
eps5s  =  Gradients.eps5s*1E-6;
eps6   =  Gradients.eps6*1E-6;
eps7   =  Gradients.eps7*1E-6;

%--------------------------------------------------------------------------
%% TIMING

% capital deltas
D2   =   t22 - t21;
D3   =   t32 - t31;
D4   =   t42 - t41;
D5rp =   TE  - t5rp;
D5s  =   TE  - t5s;
D71  =   TE  - t71;
D75  =   D71;

% tau-s
tau11       =   d1^3;
tau22       =   d2^2*(D2-d2/3)+eps2^3/30-d2*eps2^2/6;
tau23       =   d2*d3*D3;
tau24       =   d2*d4*D4;
tau25s      =   d5s*d2*D2/2;
tau25rp     =   d5rp*d2*D2/2;
tau33       =   d3^2*(D3-d3/3)+eps3^3/30-d3*eps3^2/6;
tau34       =   d3*d4*D4;
tau35s      =   d5s*d3*D3/2;
tau35rp     =   d5rp*d3*D3/2;
tau44       =   d4^2*(D4-d4/3)+eps4^3/30-1/6*d4*eps4^2;
tau45s      =   d5s*d4*D4/2;
tau45rp     =   d5rp*d4*D4/2;
tau55rp     =   d5rp^2*(D5rp-d5rp/3)+eps5rp^3/30-d5rp*eps5rp^2/6;
tau55s      =   d5s^2*(D5s-d5s/3)+eps5s^3/30-d5s*eps5s^2/6;
tau5rp71    =   d5rp*(d7*(D75-d7/4)+eps7^2/12-d7*eps7/2);
tau5s71     =   d5s*(d7*(D75-d7/4)+eps7^2/12-d7*eps7/2);
tau6m71     =   eps6/4*(d7*D71-eps7^2/60);
tau7171     =   1/4*(d7^2*(D71-d7/3)+eps7^3/30-d7^2*eps7/2);
tau7m7mplus =   (res/2-1)*(d7^3/12+eps7/60+d7^2*eps7/4-d7*eps7^2/12); 

% taus with summation
tau5rp6m = 0;
tau5s6m  = 0;
tau6m6m  = 0;

for m=1:res/2  
    
    t6i =   Gradients.t6(m);
    D6m =   TE  - t6i;
    
    tau5rp6m = tau5rp6m + eps5rp*d5rp*(D6m-eps6);
    tau5s6m  = tau5s6m  - eps5s*d5s*(D6m-eps6);
    tau6m6m  = tau6m6m  + eps6^2 * ((2*m-1)*D6m-(67*m/30-1)*eps6);
    
end

%--------------------------------------------------------------------------
%% b-matrix terms

% diagonal terms
brr = gsq*(Gdr^2*tau22+2*Gdr*Gcr*tau23+Gcr^2*tau33+Grdp^2*tau55rp+...
    2*Grdp*Gro*tau5rp71+Gro^2*(tau7171+tau7m7mplus));

bpp = gsq*(Gdp^2*tau22+2*Gdp*Gcp*tau23+Gcp^2*tau33+Gpdp^2*tau55rp+...
    2*Gpdp*Gpe*tau5rp6m+Gpe^2*tau6m6m);

bss = gsq*(14/3*Gsl^2*tau11+Gds^2*tau22+2*Gds*Gcs*tau23+Gds*Gsl2*tau24+...
    Gds*Grf*tau25s+Gcs^2*tau33+Gcs*Gsl2*tau34+2*Gcs*Grf*tau35s+...
    Gsl2*tau44/4+Gsl2*Grf*tau45s+Grf^2*tau55s/4);

% off-diagonal terms
brp = gsq*(Gdr*Gdp*tau22+(Gdr*Gcp+Gcr*Gdp)*tau23+Gcr*Gcp*tau33+Grdp*Gpdp*tau55rp+...
    +Grdp*Gpe*tau5rp6m+Gro*Gpdp*tau5rp71+Gpe*Gro*tau6m71);
bpr = brp;

brs = gsq*((Gdr*Gcs+Gcr*Gds)*tau23+1/2*Gdr*Gsl2*tau24+1/2*Gdr*Grf*tau25s+...
    Gcr*Gcs*tau33+Gcr*1/2*Gsl2*tau34+1/2*Gcr*Grf*tau35s+1/2*Grdp*Grf*tau55rp+...
    1/2*Gro*Grf*tau5s71);
bsr = brs;

bsp = gsq*(Gdp*Gds*tau22+(Gds*Gcp+Gdp*Gcs)*tau23+1/2*Gdp*Gsl2*tau24+...
    1/2*Gdp*Grf*tau25s+Gcp*Gcs*tau33+1/2*Gcp*Gsl2*tau34+1/2*Gcp*Grf*tau35s+...
    1/2*Gpdp*Grf*tau55rp+1/2*Gpe*Grf*tau5s6m);
bps = bsp;

% matrix
b = [brr, brp, brs;...
     bpr, bpp, bps;...
     bsr, bsp, bss];
 
end
