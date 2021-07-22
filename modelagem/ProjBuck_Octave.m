%% PREPARACAO
clear all;          % limpeza das variaveis
close all;          % fecha janelas abertas
clc;                % limpeza do console
pkg load control;

set(0, 'defaultlinelinewidth', 1.5);

xline = @(xval, varargin) ...
         line([xval xval], ylim, varargin{:});

function retval = manual_c2d(G, tsamp, method='f')
   z = tf('z', tsamp);
   
   % aproximação forward/backward/trapezoidal
   if method=='f'       s = (z-1)/tsamp;
   elseif method=='b'   s = 1/tsamp * (z-1)/z;
   else                 s = 2/tsamp * (z-1)/(z+1);
   end

   numz = 0;
   denz = 0;
   [num den] = tfdata(G,'v');

   for k=1:length(num)
      numz = numz*s + num(k);
   end
   for k=1:length(den)
      denz = denz*s + den(k);
   end

   retval = minreal(numz/denz);
end

function discrete_rec_eq(Gz, resolution=3, input='e', output='u')
   [a b] = tfdata(Gz,'v'); % coeficientes
   n = length(a)-1;        % ordem do numerador
   m = length(b)-1;        % ordem do denominador
   
   recstr = sprintf(['\n' output '(k-0) = ']);
   numstr = ['%+.' int2str(resolution) 'e'];
   
   format = ['...\n   ' numstr ' * ' input '(k-%i)'];
   for k=0:n
      if a(k+1) != 0
         recstr = [recstr sprintf(format, a(k+1)/b(1), k)];
      end
   end
   
   format = ['...\n   ' numstr ' * ' output '(k-%i)'];
   for k=1:m
      if b(k+1) != 0
         recstr = [recstr sprintf(format, -b(k+1)/b(1), k)];
      end
   end
   
   printf('%s;\n', recstr);
end

function retval = sat(value, min_v, max_v)
   % bloco de saturação
   retval = min(max(value, min_v), max_v);
end

f_r = 50;    % frequencia do sinal de entrada
f_s = 15000; % frequencia de chaveamento

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CIRCUITOS ANALÓGICOS

f_cb = sqrt(f_r*f_s)
w_cb = 2*pi*f_cb;           % freq. natural
z_cb = 0.46;                % amortecimento
num = [0 0 w_cb^2];
den =[1 2*w_cb*z_cb w_cb^2];
G = tf(num, den)            % FT do Buck

alpha = 0.5549;
f_sk = alpha*sqrt(f_s*f_cb)
w_sk = 2*pi*f_sk;           % freq. natural
z_sk = 0.4;                 % amortecimento
num = [0 0 w_sk^2];
den = [1 2*w_sk*z_sk w_sk^2];
H = tf(num, den)            % FT do Sallen-Key

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FILTRO DIGITAL

f_dg = alpha*sqrt(f_s*f_sk)
w_dg = 2*pi*f_dg;           % freq. natural
z_dg = 0.2;                 % amortecimento
num = [0 0 w_dg^2];
den = [1 2*w_dg*z_dg w_dg^2];
F = tf(num, den)            % FT do filtro digital

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RESPOSTA EM FREQUENCIA

W = logspace(1, 6, 100);
[Mg, Pg, _] = bode(G,W);
[Mh, Ph, _] = bode(H,W);
[Mf, Pf, _] = bode(F,W);
W = W/(2*pi);

fig = figure();
hold on;

subplot(2,1,[1]);
semilogx(W, 20*log10(1/Mg), ...
         W, 20*log10(1/Mh), ...
         W, 20*log10(1/Mf));
xlim([4e1 2e4]);
ylim([-120 -20]);
xline(f_r, 'linestyle', '-.', 'color', [0 0.7 0]);
xline(f_s, 'linestyle', '-.', 'color', [.7 0 .7]);
ylabel('Magnitude [dB]');
set(gca,'fontsize',12);
grid on;

subplot(2,1,2);
semilogx(W, Pg, W, Ph, W, Pf);
xlim([4e1 2e4]);
ylim([-200 20]);
xline(f_r, 'linestyle', '-.', 'color', [0 0.7 0]);
xline(f_s, 'linestyle', '-.', 'color', [.7 0 .7]);
ylabel('Fase [deg]');
xlabel('Frequência [Hz]');
lbl = legend('G(s)','H(s)','F(s)', 'fr','fs');
set(lbl,'fontsize',14);
set(gca,'fontsize',12);
grid on;

print(fig, 'diag_bode.png', '-r300');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RESPOSTA AO DEGRAU

tdelta = 1e-5;    % passo de tempo
tfinal = 2.5e-3;  % tempo final
[xG, t] = step(G, tfinal, tdelta);
[xH, _] = step(H, tfinal, tdelta);
[xF, _] = step(F, tfinal, tdelta);
t = t*1000; % mudança de unidade

fig = figure();
hold on;
grid on;
plot(t, xG);
plot(t, xH);
plot(t, xF);
ylabel('Tensao [V]');
xlabel('Tempo [ms]');
lbl = legend('G(s)','H(s)', 'F(s)');
set(lbl,'fontsize',14);
set(gca,'fontsize',12);

print(fig, 'resp_degrau.png', '-r300');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PROJETO PID

% Controlador tipo PID
%    C_pid(s) = (Kp + Ki/s + Kd s)
%    C_pid(s) = (Kd s^2 + Kp*s + Ki) / s

% Critério de estabilidade de Routh-Hurwitz
%   R(s) -->(+)-->[ K*G(s) ]--+--> C(s)
%            ^                |
%            |                |
%            +----[ -H(s) ]<--+
%
%     T(s) = C(s)/R(s)
%     T(s) = K G(s) / [1 + K G(s) H(s)]
%     T(s) = [K NG(s)/DG(s)] / [1 + [K NG(s) NH(s)]/[DG(s) DH(s)]]
%     T(s) = [K NG(s) DH(s)] / [DG(s) DH(s) + K NG(s) NH(s)]
%  tomando T(s) = NT(s)/DT(s), aplica-se o critério em
%     DT(s) = DG(s) DH(s) + K NG(s) NH(s)


[numG denG] = tfdata(G,'v');
[numH denH] = tfdata(H,'v');
denT_p1 = tf(denG,1) * tf(denH,1); % termo s/ a constante K
denT_p2 = tf(numG,1) * tf(numH,1); % termo c/ a constante K
%---------------------------------------------------------
%  denT(s) = [ 1 1.506e4 2.378e8 1.088e12 4.676e15*(1+K) ]
%---------------------------------------------------------
%  s^4  |     1      2.378e8    4.676e15(1+K)  |
%  s^3  |  1.506e4   1.088e12         0        |
%---------------------------------------------------------
%  s^2  |  1.656e8   4.676e15(1+K)    0        | K > -1
%  s^1  | (6.628e11 -      0          0        | K < 1.5588
%       |  K*4.252e11)                         |
%  s^0  |  4.676e15(1+K)   0                   |
%---------------------------------------------------------
Ku = 1.5588
%step(feedback(Ku*G, H), tfinal, delta);
%xticks(0:2e-4:2*tfinal);
Tu = 0.00095 - 0.0002

Kp = 0.3*Ku
Ki = 1.5*Ku/Tu
Kd = 0.05*Ku*Tu

C = tf([Kd Kp Ki], [1 0])
Gcmf = feedback(C*G, H*F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DISCRETIZAÇÃO DA PLANTA

VCC = 3.3;
VSOURCE = 24;
Ka = VSOURCE;     % ganho do acionamento
Kt = VCC/VSOURCE; % ganho do transdutor

NC = 2;        % número de ciclos
T  = 1/f_r;    % período da referência
N  = f_s/f_r;  % amostras por período
TT = NC*T;     % duração da simulação
NT = NC*N;     % amostras na simulação

NP = 50;      % resolução da superdiscretização
dn = T/N;      % intervalo entre amostras
dt = dn/NP;    % intervalo de superdiscretização
n = 0:dn:(NC*T-dn);  % instantes de amostragem
t = 0:dt:(NC*T-dt);  % tempo superdiscretizado

% discretização e superdiscretização
Cz = manual_c2d(C, dn);       % controlador digital
Fz = manual_c2d(F, dn,'t');   % filtro digital
Gz = manual_c2d(G*Ka, dt);    % conversor buck
Hz = manual_c2d(H*Kt, dt);    % filtro analógico
PIDz = [
   manual_c2d(tf(Kp), dn)     % parcela proporcional
   manual_c2d(Ki/tf('s'), dn) % parcela integral
   manual_c2d(Kd*tf('s'), dn) % parcela derivativa
];

% equações recursivas
prec = 3; % número de casas decimais (precisão)
discrete_rec_eq(PIDz(1,:), prec, 'e', 'up'); % proporcional
discrete_rec_eq(PIDz(2,:), prec, 'e', 'ui'); % integrador
discrete_rec_eq(PIDz(3,:), prec, 'e', 'ud'); % derivativo
discrete_rec_eq(Fz, 2*prec, 'yq', 'yf');  % filtro digital
discrete_rec_eq(Gz, 2*prec, 'dc', 'y');  % conversor buck
discrete_rec_eq(Hz, 2*prec, 'y', 'ym'); % filtro analógico

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIMULAÇÃO EM TEMPO DISCRETO

NBITS = 10;
RES_ADC = 2^NBITS;

% variáveis discretas
x  = zeros(NT,1);   % referência
e  = zeros(NT,1);   % erro
u  = zeros(NT,1);   % ação de controle
up = zeros(NT,1);   % ação proporcional
ui = zeros(NT,1);   % ação integradora
ud = zeros(NT,1);   % ação derivativa
yq = zeros(NT,1);   % medição quantizada
yf = zeros(NT,1);   % medição filtrada
dc = zeros(NT,1);   % duty cicle

% variáveis pseudo-contínuas
y  = zeros(NT*NP,1); % tensão de saída
ym = zeros(NT*NP,1); % saída medida

% geração do sinal de referência
for k=1:NT
   idx = mod(k, N);
   if idx < (9e-3/dn)
      x(k) = 6;
   elseif idx < (10e-3/dn)
      x(k) = -102 + (12/1e-3)*(idx*dn);
   elseif idx < (19e-3/dn)
      x(k) = 18;
   else
      x(k) = +246 - (12/1e-3)*(idx*dn);
   end
end

%plot(n,x);
%ylim([0 24]);
%grid on;

% construção das curvas
for k=3:NT
   % Etapa de amostragem (ADC -> SNH) e quantização
   adc_in = round(ym((k-1)*NP) * RES_ADC);
   yq(k) = adc_in / (Kt*RES_ADC);
   % Filtro digital e erro
   yf(k) = ...
      +2.441204e-01 * yq(k-0)...
      +9.764818e-01 * yq(k-1)...
      +1.464723e+00 * yq(k-2)...
      +9.764818e-01 * yq(max(k-3,1))...
      +2.441204e-01 * yq(max(k-4,1))...
      -1.283287e+00 * yf(k-1)...
      -2.597685e-01 * yf(k-2)...
      -6.696768e-01 * yf(max(k-3,1))...
      -6.931951e-01 * yf(max(k-4,1));
   e(k) = x(k) - yf(k);
   % Ações de controle individuais
   up(k) = +4.676e-01 * e(k);
   ui(k) = +2.078e-01 * e(k) + ui(k-1);
   ud(k) = +8.768e-01 * (e(k) - e(k-1));
   % Ação de controle conjunta com saturação
   ui(k) = sat(ui(k), -VSOURCE, +VSOURCE);
   u(k) = sat(up(k) + ui(k) + ud(k), 0, VSOURCE);
   % Duty cicle
   dc(k) = u(k)/Ka;
   % Resposta do sistema (simulação da planta real)
   for j=1:NP
      kc = (k-1)*NP + j;
      % Resposta da planta
      y(kc) = ...
         +1.263309e-03 * dc(k-0)...
         +1.993325e+00 * y(kc-1)...
         -9.933779e-01 * y(kc-2);
      % Saída do Sallen-Key
      ym(kc) = ...
         +3.860035e-05 * y(kc-0)...
         +1.986596e+00 * ym(kc-1)...
         -9.868767e-01 * ym(kc-2);
      % Saturação do Amp-Op
      ym(kc) = sat(ym(kc),0,VCC);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXIBIÇÃO GRÁFICA

figure();
grid on;
hold on;
stairs(n,x);
plot(t,y);
xlabel('Tempo [s]');
ylabel('Tensão [V]');
lbl = legend('x[n]', 'y(t)', ...
      'location', 'southeast');
set(lbl,'fontsize',14);
set(gca,'fontsize',12);
print('simu_x_y.png', '-r300');
close;


figure();
grid on;
hold on;
stairs(n,e);
stairs(n,u);
xlabel('Tempo [s]');
ylabel('Tensão [V]');
lbl = legend('e[n]', 'u[n]', ...
      'location', 'southeast');
set(lbl,'fontsize',14);
set(gca,'fontsize',12);
print('simu_e_u.png', '-r300');
close;


figure();
grid on;
hold on;
stairs(n,up);
stairs(n,ui);
stairs(n,ud);
xlabel('Tempo [s]');
ylabel('Tensão [V]');
lbl = legend('up[n]','ui[n]','ud[n]', ...
      'location', 'northwest');
set(lbl,'fontsize',14);
set(gca,'fontsize',12);
print('simu_u_pid.png', '-r300');
close;


figure();
grid on;
hold on;
stairs(n,yq);
stairs(n,yf);
plot(t,y);
xlim([1.4*T 1.6*T]);
ylim([5 20]);
xlabel('Tempo [s]');
ylabel('Tensão [V]');
lbl = legend('yq[n]', 'yf[n]','y(t)', ...
      'location', 'southeast');
set(lbl,'fontsize',14);
set(gca,'fontsize',12);
print('simu_q_f.png', '-r300');
close;