%% by Bellila Ahmed Nassim, student at ENST
%On exchange at Université Laval
%Course UL GEL-7029 Prictive control under supervision of Prof. André Desbiens

%see Readme.txt to understand more the process

%%
clc;clear all;

%% Sim parameters
% sampling period
T=0.5;
% number of points to simulate
Npts=50;
 
% plant
numdp=[0 0.4]; % WITHOUT the deadtime from the ZOH
dendp=[1 -0.8];
 
% matrix to record the simulation data 
simdata=zeros(Npts,3);
% initialization
xp=zeros(max(length(dendp),length(numdp))-1,1); % plant states
% DMC
g=[ 0; 0.4; 0.7; 0.9; 1.2; 1.3; 1.5; 1.6; 1.66; 1.73; 1.78; 1.8; ...
 1.9; 2.0; 2.1; 2.0; 2.1; 1.9; 2.0; 2.1; 2.1; 1.9; 2.0; 2.0; 1.9];
 % g = recorded unit plant step response
 % g(1) = g(k=1); etc.
Lambda=0; % weight on the control moves
Hp=10; % prediction horizon
Hc=1; % control horizon 
% set-point for the whole simulation
cons=[0; 1*ones(Npts-1,1)];

%%
%Here we go 1 variable initialisation
%G*********************
G = zeros(Hp, Hc);
for j=1:Hp
    G(j,1) = g(j);
end

if Hc>1
    for i=2:Hc
        G(:,i) = [zeros(i-1,1); g(1:Hp-i+1)];
    end    
end
K = inv(G'*G + Lambda*eye(length(G'*G)))*G';
K = K(1,:);

N = 15;
f0=zeros(Hp,N-1);


%% matrix f0*************
for j=1:Hp
    for i=N:-1:2
        f0(j,N-i+1) = g(i+j-1) - g(i-1);
    end
end


%% initialisation        
Vupast = zeros(N-1,1);
Vupast1 = zeros(N-1,1);
Deltaupast = Vupast - Vupast1;
upast = 0;
%vecteur de valeurs de u (N past values and simulation values)
Vu = zeros(Npts+15,1);


%% simulation loop
for i=16:1:Npts+15
    % plant output y
    [y,xp]=filter(numdp,dendp,upast,xp);
    %%
    %here we go 2 calculate the controle u(k)
    for l=N-1:-1:1
        Vupast(N-l) = Vu(i-l);
        Vupast1(N-l) = Vu(i-l-1);
        Deltaupast = Vupast - Vupast1;
    end
    %f = y(k) + f0*Deltau
    f = f0*Deltaupast+y;
    u = K*(cons(i-15)-f) + upast;   
    upast = u;
    Vu(i)=u;
    % data recording
    simdata(i,1:3)=[cons(i-15) y u];
end


%% plot data
%Pulling back the N past values of vec simdata
simdata = simdata(16:Npts+15,:);

%figures
figure(1);
t=T*(0:Npts-1)';
subplot(211);
stairs(t,simdata(:,3),'k-');
ylabel('u');
subplot(212);
plot(t,cons,'k--',t,simdata(:,2),'r-');
ylabel('r (black) and y (red)');
xlabel('Time [sec]'); 