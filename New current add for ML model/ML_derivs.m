function dydt = ML_derive(t, y, I, coefs)
% The derivs function

% Uncomment the section for simulation
%% A. Original model parameter.

% persistent V;
% persistent n;
% persistent Ser;
% persistent Sef;
% persistent p;

% V = y(1,:);
% n = y(2,:);
% Ser = y(3,:);
% Sef = y(4,:);

% dvdt=1.0./coefs.C*(I-coefs.gL.*(V-coefs.EL) - ...
%     coefs.gNa.*(1.0./(1.0 + exp((coefs.Vhalfm - V)/coefs.km ))).*(V-coefs.ENa) - ...
%     coefs.gK.*n.*(V-coefs.EK) + ...
%     coefs.gEPSP.*(Sef-Ser).*(coefs.E_EPSP-V));

% dndt=(1.0./(1.0+exp((coefs.Vhalfn - V)/coefs.kn ))-n)./exp(-0.07*V-3);

% dsdt=-Ser/coefs.tauEPSPr;
 
% dfdt=-Sef/coefs.tauEPSPf;

% dydt=[dvdt;dndt;dsdt;dfdt];

%% B. Fast persistent Na+ current testing parameter, test whether model can firing when it appear in model.

% persistent V;
% persistent n;
% persistent p;

% V = y(1,:);
% n = y(2,:);
% p = y(3,:);

% dvdt=1.0./coefs.C*(I-coefs.gL.*(V-coefs.EL) - ...
%     coefs.gNa.*(1.0./(1.0 + exp((coefs.Vhalfm - V)/coefs.km ))).*(V-coefs.ENa) - ...
%     coefs.gK.*n.*(V-coefs.EK) + ...
%     coefs.gNaP.*p.*(coefs.ENa-V));

% dndt=(1.0./(1.0+exp((coefs.Vhalfn - V)/coefs.kn ))-n)./exp(-0.07*V-3);

% dpdt = (1 ./ (1 + exp((-54 - V) ./ 9)) - p) ./ 0.8;

%dydt=[dvdt;dndt;dpdt];

%% B_1. Fast persistent Na+ current simulation for model.

% persistent V;
% persistent n;
% persistent Ser;
% persistent Sef;
% persistent p;

% V = y(1,:);
% n = y(2,:);
% Ser = y(3,:);
% Sef = y(4,:);
% p = y(5,:);


% dvdt=1.0./coefs.C*(I-coefs.gL.*(V-coefs.EL) - ...
%     coefs.gNa.*(1.0./(1.0 + exp((coefs.Vhalfm - V)/coefs.km ))).*(V-coefs.ENa) - ...
%     coefs.gK.*n.*(V-coefs.EK) - ...
%     coefs.gNaP.*p.*(V-coefs.ENa)+...
%     coefs.gEPSP.*(Sef-Ser).*(coefs.E_EPSP-V));

% dndt=(1.0./(1.0+exp((coefs.Vhalfn - V)/coefs.kn ))-n)./exp(-0.07*V-3);

% dsdt=-Ser/coefs.tauEPSPr;
 
% dfdt=-Sef/coefs.tauEPSPf;

% dpdt = (1 ./ (1 + exp((-54 - V) ./ 9)) - p) ./ 0.8;

% dydt=[dvdt;dndt;dsdt;dfdt;dpdt];

%% C.M-current testing parameter, test whether model can firing when it appear in model.

% persistent V;
% persistent n;
% persistent w; 

% V = y(1,:);
% n = y(2,:);
% w=y(3,:);

% dvdt=1.0./coefs.C*(I-coefs.gL.*(V-coefs.EL) - ...
%     coefs.gNa.*(1.0./(1.0 + exp((coefs.Vhalfm - V)/coefs.km ))).*(V-coefs.ENa) - ...
%     coefs.gK.*n.*(V-coefs.EK)) - ...
%     coefs.Mcur.* w.*(V-coefs.EK);

% dndt=(1.0./(1.0+exp((coefs.Vhalfn - V)/coefs.kn ))-n)./exp(-0.07*V-3);

% dwdt = (1./(1+exp(-(V+35)./10)) - w)./(400./(3.3*exp((V+35)./20))+exp(-(V+35)./20));

% dydt=[dvdt;dndt;dwdt];

%% C_1. M-current simulation for model.

% persistent V;
% persistent n;
% persistent Ser;
% persistent Sef;
% persistent w; 

% V = y(1,:);
% n = y(2,:);
% Ser = y(3,:);
% Sef = y(4,:);
% w=y(5,:);

% dvdt=1.0./coefs.C*(I-coefs.gL.*(V-coefs.EL) - ...
%     coefs.gNa.*(1.0./(1.0 + exp((coefs.Vhalfm - V)/coefs.km ))).*(V-coefs.ENa) - ...
%     coefs.gK.*n.*(V-coefs.EK) - ...
%     coefs.Mcur.* w.*(V-coefs.EK)+...
%     coefs.gEPSP.*(Sef-Ser).*(coefs.E_EPSP-V));

% dndt=(1.0./(1.0+exp((coefs.Vhalfn - V)/coefs.kn ))-n)./exp(-0.07*V-3);

% dsdt=-Ser/coefs.tauEPSPr;

% dfdt=-Sef/coefs.tauEPSPf;

% dwdt = (1./(1+exp(-(V+35)./10)) - w)./(400./(3.3*exp((V+35)./20))+exp(-(V+35)./20));

% dydt = [dvdt;dndt;dsdt;dfdt;dwdt];

%% Optional: Testing difference with original model and model with additional current.

% persistent V;
% persistent n;

% V = y(1,:);
% n = y(2,:);

% dvdt=1.0./coefs.C*(I-coefs.gL.*(V-coefs.EL) - ...
%     coefs.gNa.*(1.0./(1.0 + exp((coefs.Vhalfm - V)/coefs.km ))).*(V-coefs.ENa) - ...
%     coefs.gK.*n.*(V-coefs.EK));
 
% dndt=(1.0./(1.0+exp((coefs.Vhalfn - V)/coefs.kn ))-n)./exp(-0.07*V-3);

% dydt=[dvdt;dndt];












