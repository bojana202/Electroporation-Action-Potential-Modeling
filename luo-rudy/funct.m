function dy=fun_LR91(t,y,flags,data)

% State variables
V = y(1); H = y(2); m = y(3); J = y(4); d = y(5); f = y(6); x = y(7); Ca_i = y(8);

% Simulus current
In = data.Is*(t<=data.ts);

% Ion channel currents
[INa,am,bm,aH,bH,aj,bj] = comp_INa(V,m,H,J,data);
[Isi,ad,bd,af,bf] = comp_Isi91(V,d,f,Ca_i,data);  
[IK,ax,bx] = comp_IK(V,x,data);
[IK1] = comp_IK1(V,data);
[IKp] = comp_IKp(V,data);

% Background current
Ib = data.Gb_max*(V+59.87);

% Electroporation current
Iep = data.Gep*(V-0);

% Compute derivatives of the state variables
dV = In-(INa+IK+Ib+IKp+IK1+Isi+Iep);
dH = aH.*(1-H)-bH.*H;
dm = am.*(1-m)-bm.*m;
dJ = aj.*(1-J)-bj.*J;
dD = ad.*(1-d)-bd.*d;
df = af.*(1-f)-bf.*f;
dx = ax.*(1-x)-bx.*x;
dCa_i = -1e-4*Isi + 0.07*(1e-4-Ca_i); % TO BE UPGRADED WIH EP

% RETURN DERIVATIVES
dy = [dV; dH; dm; dJ; dD; df; dx; dCa_i];

function [INa,am,bm,ah,bh,aj,bj] = comp_INa(V,m,H,J,data)
% INa ... the fast sodium current in mammalian ventricular cells

ENa = data.RTF*log(data.Na_o./data.Na_i); % Nernst potential of Na, mV
INa = data.GNa_max*m.*m.*m.*H.*J.*(V-ENa);

a=1-1./(1+exp(-(V+40)/0.24)); % sigmoidal function: a = 1 for V < -40, a = 0 for V > -40
ah= a.*0.135.*exp((80+V)./(-6.8));
bh= (1-a)./(0.13*(1+exp((V+10.66)/(-11.1)))) +(a).*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));

aj =  a.*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V)).*(V+37.78)./(1+exp(0.311*(V+79.23)));
bj= (1-a).*(0.3*exp(-2.535e-7*V)./(1+exp(-0.1*(V+32))))+(a).*(0.1212*exp(-0.01052*V)./(1+exp(-0.1378*(V+40.14))));

am = 0.32*(V+47.13)/(1-exp(-0.1*(V+47.13)));
bm = 0.08*exp(-V/11);

function [Is,ad,bd,af,bf] = comp_Isi91(V,d,f,Ca_i,data)
% Isi ... the slow inward current
 
Esi = 7.7 - 13.0287*log(Ca_i./data.Ca_o); % In the L-R paper there is no normalization on Ca_o
Is =  data.Gsi_max*d.*f.*(V-Esi);

ad = 0.095*exp(-0.01*(V-5))./(1+exp(-0.072*(V-5)));
bd = 0.07*exp(-0.017*(V+44))./(1+exp(0.05*(V+44)));

af = 0.012*exp(-0.008*(V+28))./(1+exp(0.15*(V+28)));
bf = 0.0065*exp(-0.02*(V+30))./(1+exp(-0.2*(V+30)));

function [IK,ax,bx] = comp_IK(V,x,data)
% IK ... Time-dependent potassium current

Xi=1*(V<=-100)+(V > -100).*(2.837*(exp(0.04*(V+77))-1)./((V+77).*exp(0.04*(V+35))));

EK =data.RTF.*log((data.K_o+data.PNa_K*data.Na_o)/(data.K_i+data.PNa_K*data.Na_i));
IK = data.GK_max.*Xi.*x.*(V-EK);

ax = 0.0005*exp(0.083*(V+50))./(1+exp(0.057*(V+50)));
bx = 0.0013*exp(-0.06*(V+20))./(1+exp(-0.04*(V+20)));

function [IK1] = comp_IK1(V,data)
% IK1 ... Time-independent potassium current

EK1 = data.RTF*log(data.K_o./data.K_i);

ak1 = 1.02./(1+exp(0.2385*(V-EK1-59.215)));
bk1 = (0.49124*exp(0.08032*(V-EK1+5.476))+exp(0.06175*(V-EK1-594.31)))./...
    (1+exp(-0.5143*(V-EK1+4.753)));

k1_inf = ak1./(ak1+bk1);
IK1 =data.GK1_max*k1_inf.*(V-EK1);

function [IKp] = comp_IKp(V,data)
% IKp ... Plateau potassium current
EKp = data.RTF*log(data.K_o./data.K_i); % Same as EK1

kp_inf = 1./(1+exp((7.488-V)/5.98));
IKp = data.GKp_max*kp_inf*(V-EKp);
