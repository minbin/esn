clear;
format long;

global E Vp Vi Vg tp ti td Rm Rg C1 C2 C3 C4 C5 k a1 Ub U0 Um a b G
global A K T eps

Vp = 3;
Vi = 11;
Vg = 10;
E  = .2;
tp = 6;
ti = 100;
td = 12;
k  = 0.5;
Rm = 209;
a1 = 6.6;
C1 = 300;
C2 = 144;
C3 = 100;
C4 = 80;
C5 = 26;
Ub = 72;
U0 = 4;
Um = 94;
Rg = 180;
a  = 7.5;
b  = 1.772;
G = 0;

A = 1000;
K = 2000;
T = 200;
eps = 1;

% delta_pulses([0 0 0 0 0 0],100,1000,20)
gaussian_pulses([0 0 0 0 0 0],eps,K,A,T,'n')

function delta_pulses(inits,N,A,T)
  % Parameters: (initial conditions, number of iterations, kick, inter-kick time)

  % pre-allocate vector for initial conditions
  I = cat(1, inits, zeros(N, size(inits,2)));
  U = zeros(N, size(inits,2));

  % store glucose and time
  Y = [];
  tt = [];

  for kk = 1:N
    tspan = [(kk-1)*T kk*T];
    y0 = I(kk,:);
    
    [t,y] = ode23s(@(t,y) Ultradian(t,y,'delta'), tspan, y0);

    % implement delta kick
    y(end,3) = y(end,3) + A;

    % store solutions
    Y = [Y; y];
    tt = [tt; t];

    % store ultradian output
    U(kk,:) = Ultradian(tspan, y0, 'delta');

    % modify initial conditions
    I(kk+1,:) = Y(end,:);
  end

  figure(1)
  plot(tt,Y(:,3),'k','LineWidth',3)
  set(gca,'fontsize',20)
  ylabel('glucose')
  xlabel('time')

  writematrix(tt, 'figDP_tt.csv')
  writematrix(Y(:,3), 'figDP_Y3.csv')
  writematrix(Y, 'figDP_Y.csv')
  writematrix(I, 'figDP_I.csv')
  writematrix(U, 'figDP_U.csv')
end

function gaussian_pulses(inits,eps,K,A,T,graph_pulse)
  % Parameters: (initial conditions, epsilon, K, inter-kick time)
  
  % final time
  TT = K*T-1;
  
  [t,y] = ode23s(@(t,y) Ultradian(t,y,'square'),[0,TT],inits);

  figure(1)
  plot(t,y(:,3),'k','LineWidth',3);
  set(gca,'fontsize',20)
  ylabel('glucose')
  xlabel('time')
  
  if strcmp(graph_pulse,'y')
    tt = 0:0.01:TT-1;
    y2 = 0;
    for n = 1:K 
      y2 = y2 + (A/(eps*sqrt(pi)))*exp(-((tt-n*T)/eps).^2);
    end

    figure(2)
    plot(tt,y2,'k','LineWidth',3)
    set(gca,'fontsize',20)
    ylabel('square pulse')
    xlabel('time')
  end
  
  writematrix(t, 'figGp_t.csv')
  writematrix(y, 'figGp_y.csv')
  writematrix(y(:,3), 'figGp_y3.csv')
end
