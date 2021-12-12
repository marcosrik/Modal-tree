%===========================================%
%------ Modal Analysis of Fractal trees ----%
%-------------------------------------------%
%        FEM frame element algorithm        % 
%-------------------------------------------%
% Autor: Marcos Antonio R. Ezequiel         %
%-------------------------------------------%
% Data: 11/12/2021                          %
%===========================================%

%% Initialization
clc; clear all; close all;

%% --------------------------Input-------------------------------------%
% Base diameter
D_base = 18e-2; % [m]
% inclination of the branches 
alpha = 20; % [°]
% Parameter that determines 
%the length of the branches through L = D^Beta
Beta = 1.5;
% Area reducing factor A1 = A0*lambda
lambda = .5;
% Density
rho = 805; % [kg/m^3]
% Modulus of Elasticity 
E = 11.3e9; %[Pa]
% Poisson coefficient
nu = 0.38;
% Number of branch levels
Levels = 4;
% Number of modes
nmodes = 10;
% Scale Factor for the plot
scale = 90;

%% --------------------------Mesh--------------------------------------%
L = zeros(1,Levels);
D = zeros(1,Levels);
D_i = D_base;

% Coordinates matrix coord[x , y]
coord = zeros(2^Levels,2);
coord(1,:) = [0,0];
coord_i = [0,0];
cc = 2;

%Incidence Matrix inci[Tmat,Tgeo,no1,no2]
nel = 2^Levels-1;
inci = zeros(nel,4);
inci(:,1) = ones(nel,1); 
e = 1;
nodes_lev = 1;

theta = 0;
for i = 1:Levels
%Properties of each level
D(i) = D_i;
L(i) = D(i)^(1/Beta);
D_i = D(i)*sqrt(lambda);
%Coordinates and incidence matrix
nnodes_lev = 2^(i-1);
for ii = 1:nnodes_lev
    coord(cc,:) = [sind(theta(ii))*L(i),...
                  cosd(theta(ii))*L(i)]+coord_i(ii,:);           
    inci(e,2:4) = [i,nodes_lev(ii),cc];
    
    cc = cc+1;
    e = e+1;
end
nodes_lev = repelem(cc-nnodes_lev:cc-1,1,2);
coord_i = coord(nodes_lev,:);
theta = repelem(theta,1,2)+repmat([-alpha,alpha],1,nnodes_lev);
end

%% --------------------------Mesh Plot---------------------------------%
nnodes = length(coord(:,1));
vnodes = [1:nnodes]';
nnel = 2;      %Nodes per element
X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
for iel=1:nel
    for i=1:nnel
        X(i,iel)=coord(inci(iel,i+2),1);
        Y(i,iel)=coord(inci(iel,i+2),2);
    end
end
plot(X,Y,'k')
title('Finite Element Mesh') ;
fill(X,Y,'w') 
axis off
axis equal
%Legend
text(coord(:,1),coord(:,2),int2str(vnodes(:,1)),'fontsize',8,'color','k');
for i = 1:nel
    text(sum(X(:,i))/nnel,sum(Y(:,i))/nnel,int2str(i),'fontsize',10,'color','r');
end
plotMalha(coord,inci,1,'k')
title([num2str(Levels) ' Level Tree'])


%% -------------------------Table of materials-------------------------%
% line 1: Modulus of elasticity
% line 2: Density
% line 3: Poisson coefficient

Tmat = [E; rho; nu];
    
%% -------------------------Table of geometry--------------------------%
% line 1: Area
% line 2: Moment of Inercia 
%

Tgeo = [pi*(D./2).^2;
        1/4*pi*(D./2).^4];
 

%% -------------------------Matrix Assembly----------------------------%

Kg = zeros(3*nnodes,3*nnodes);
Mg = zeros(3*nnodes,3*nnodes);

for i = 1:nel
   % Element Matrix
   node1 = inci(i,3);
   node2 = inci(i,4);
   x1 = coord(node1,1);
   y1 = coord(node1,2);
   x2 = coord(node2,1);
   y2 = coord(node2,2);
   le = sqrt((x2-x1)^2+(y2-y1)^2);
   c = (x2-x1)/le;
   s = (y2-y1)/le;
   E = Tmat(1,inci(i,1));
   rho = Tmat(2,inci(i,1));
   nu = Tmat(3,inci(i,1));
   A = Tgeo(1,inci(i,2));
   I = Tgeo(2,inci(i,2));
   % Rotational Matrix
   T = [c s 0; -s c 0; 0 0 1];
   T = [T zeros(3); zeros(3) T];
   % Stifness Matrix
   Kl = [A*E/le  0            0          -A*E/le  0            0;
         0       12*E*I/le^3  6*E*I/le^2  0      -12*E*I/le^3  6*E*I/le^2;
         0       6*E*I/le^2   4*E*I/le    0      -6*E*I/le^2   2*E*I/le
        -A*E/le  0           0            A*E/le  0            0;
         0      -12*E*I/le^3 -6*E*I/le^2  0       12*E*I/le^3 -6*E*I/le^2;
         0       6*E*I/le^2   2*E*I/le    0      -6*E*I/le^2   4*E*I/le];
   % Mass Matrix
   Ml = A*rho*le*...
        [1/3     0            0           1/6     0            0;
         0       13/35        11*le/210   0       9/70        -13*le/420;
         0       11*le/210    le^2/105    0       13*le/420   -le^2/140;
         1/6     0            0           1/3     0            0;
         0       9/70         13*le/420   0       13/35       -11*le/210;
         0      -13*le/420   -le^2/140    0      -11*le/210    le^2/105];
   
   Ke = T'*Kl*T;
   Me = T'*Ml*T;
   
   dofe = [3*node1-2, 3*node1-1,3*node1, 3*node2-2, 3*node2-1,3*node2];
   Kg(dofe,dofe) = Kg(dofe,dofe)+Ke;
   Mg(dofe,dofe) = Mg(dofe,dofe)+Me;
end


%% --------------------------Restrictions------------------------------%

FixedDofs = [3*1-2 3*1-1 3*1];
AllDofs = 1:3*nnodes;
FreeDofs = setdiff(AllDofs,FixedDofs);

%% --------------------------Solver------------------------------------%

% Modal Analysis

[Fi,D] = eig(Kg(FreeDofs,FreeDofs),Mg(FreeDofs,FreeDofs));

Freq = sqrt(diag(D))/(2*pi); % Natural Frequency [Hz]

[FN,I] = sort(Freq);

%% ----------------------Post Processing-------------------------------%
figure('Position',[100 100 1200 800])
Mode_uv_total =  zeros(nnodes,2);
for j = 1:nmodes
    jj = I(j);
    Mode = zeros(3*nnodes,1);
    Mode(FreeDofs,1) = Fi(:,jj);
    subplot(ceil(nmodes/2),2,j)
    Mode_uv = [Mode(1:3:end-2), Mode(2:3:end-1)];
    Magnitude = zeros(nel,1);
    for i = 1:nel
        node2 = inci(i,4);
        Magnitude(i) = sqrt(Mode_uv(node2,1)^2 + Mode_uv(node2,2)^2);
    end
    maxi = max(max(Magnitude));
    Magnitude = Magnitude/maxi;
    Mode_uv_total = Mode_uv_total+Mode_uv/maxi;
    coord_d = coord +scale*1e-3*Mode_uv/maxi;
    
    %Plotting undeformed state
    plot(X,Y,'k')
    axis off
    axis equal
    hold on
    
    %Plotting Results
    num = 64;
    mat = jet(num);
    colormap(mat);
    idx = round(interp1([min(Magnitude),max(Magnitude)],[1,num],Magnitude));
    Xd = zeros(nnel,nel) ;
    Yd = zeros(nnel,nel) ;
    for iel=1:nel
        for i=1:nnel
            Xd(i,iel)=coord_d(inci(iel,i+2),1);
            Yd(i,iel)=coord_d(inci(iel,i+2),2);
        end
    end    
    for iel =1:nel
        plot(Xd(:,iel),Yd(:,iel),'color',mat(idx(iel),:),'LineWidth',2)
    end
    cb = colorbar; a = get(cb); colorTitleHandle = get(cb,'Title');
    titleString = {'Relative','Displacement'};
    set(colorTitleHandle ,'String',titleString);
    a =  a.Position; %gets the positon and size of the color bar
    dx = .06; dy = 0;
    h = .1; w = .01;
    set(cb,'Position',[a(1)+dx a(2)+dy w h])% To change size
    title(['f_n = ' num2str(FN(j)) ' [Hz]'])
end

%% ----------------------Plotting all modes----------------------------%
figure('Position',[100 100 1200 800])
coord_d = coord +scale*1e-3*Mode_uv_total;
Magnitude = zeros(nel,1);
for i = 1:nel
 node2 = inci(i,4);
 Magnitude(i) = sqrt((Mode_uv_total(node2,1))^2 + (Mode_uv_total(node2,2))^2);
end
maxi = max(max(Magnitude));
Magnitude = Magnitude/maxi;
coord_d = coord +scale*1e-3*Mode_uv_total/maxi;
%Plotting undeformed state
plot(X,Y,'k')
axis off
axis equal
hold on
%Plotting Results
colormap(mat);
idx = round(interp1([min(Magnitude),max(Magnitude)],[1,num],Magnitude));
Xd = zeros(nnel,nel) ;
Yd = zeros(nnel,nel) ;
for iel=1:nel
    for i=1:nnel
        Xd(i,iel)=coord_d(inci(iel,i+2),1);
        Yd(i,iel)=coord_d(inci(iel,i+2),2);
    end
end
for iel =1:nel
    plot(Xd(:,iel),Yd(:,iel),'color',mat(idx(iel),:),'LineWidth',2)
end
cb = colorbar; a = get(cb); colorTitleHandle = get(cb,'Title');
titleString = {'Relative','Displacement'};
set(colorTitleHandle ,'String',titleString);
a =  a.Position; %gets the positon and size of the color bar
dx = .06; dy = 0;
h = .7; w = .01;
set(cb,'Position',[a(1)+dx a(2)+dy w h])% To change size
title('All modes')