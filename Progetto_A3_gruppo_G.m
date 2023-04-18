%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRUPPO G PROGETTO CONTROLLI AUTOMATICI T                               %
%                                                                        %
% TRACCIA A3                                                             %
%                                                                        %
% Giacomo Caroli         0000922118                                      %
% Lucia Gasperini        0000921439                                      %
% Lorenzo Galfano        0000934647                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc

omega_plot_min = 1e-2;
omega_plot_max = 1e6;


%% Parametri
beta = 0.1;
alfa = 0.2;
rs=1.6;
rr=1.6;
ms=0.75;
mr=0.1;
gamma=0.05;
K=900;
nse=300;
nre=600;
ue=0;


% valori di prova per il secondo punto opzionale
epsilon1=0.5;
epsilon2=0.5;

% errore a regime nullo
e_star = 0;

% attenuazione disturbo sull'uscita
A_d = 35;
omega_d_MAX = 0.15;

% attenuazione disturbo di misura
A_n = 95;
omega_n_min = 3.5*1e3;
omega_n_MAX = 3*1e6;

% Sovraelongazione massima e tempo d'assestamento all'5%
S_100_spec = 0.05;
T_a5_spec = 0.6;

%% Creazione sistema

% funzione di trasferimento
s = tf('s');

A = [rs*(1-(2*nse+nre)/K)-ms*ue-beta-alfa*ue, -((rs*nse)/K)+gamma; -((rr*nre)/K)+beta+alfa*ue, rr*(1-(nse+2*nre)/K)-mr*ue-gamma];
B = [-ms*nse-alfa*nse; -mr*nre+alfa*nse];
C = [1 0];
D = 0;

% matrice identità
I=[1 0; 0 1];

GG = zpk(C*((s*I - A)^-1)*B + D);
display(GG);

%% Diagramma di Bode

figure(1);
bode(GG,{omega_plot_min,omega_plot_max});
grid on, zoom on;



%% Regolatore statico - proporzionale senza poli nell'origine

% consideriamo il regolatore statico negativo per compensare 
% il guadagno negativo e rispettare il criterio di bode
mu = -0.01;
RR_s = mu/s;

display(RR_s);

% aggiunta di un polo ad alte frequenze per rispettare le specifiche
polo = 1 / (1+0.012*s);

% Sistema esteso
GG_e = RR_s*GG*polo;

display(GG_e);
figure(2); %grafico GG_e
bode(GG_e,{omega_plot_min,omega_plot_max});
grid on, zoom on;


%% Diagrammi di Bode di Ge con specifiche

figure(3);
hold on;

% Calcolo specifiche S% => Margine di fase
S_100= 0.05;
xi= sqrt((log(S_100))^2/(pi^2+(log(S_100))^2)); %vale circa 0,69 (smorzamento)
display(xi)
Mf_spec = xi*100; % Mf_spec = 69
display(Mf_spec)


% Specifiche su d
omega_d_min = 0.0001; % lower bound per il plot
Bnd_d_x = [omega_d_min; omega_d_MAX; omega_d_MAX; omega_d_min];
Bnd_d_y = [A_d; A_d; -150; -150];
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche su n
Bnd_n_x = [omega_n_min; omega_n_MAX; omega_n_MAX; omega_n_min];
Bnd_n_y = [-A_n; -A_n; 200; 200];
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);

% Specifiche tempo d'assestamento (minima pulsazione di taglio)
omega_Ta_min = 1e-4; % lower bound per il plot
omega_Ta_MAX = 300/(Mf_spec*T_a5_spec); % omega_c >= 300/(Mf*T^*) (tempo di assestamento del 5%)
display(omega_Ta_MAX)
Bnd_Ta_x = [omega_Ta_min; omega_Ta_MAX; omega_Ta_MAX; omega_Ta_min];
Bnd_Ta_y = [0; 0; -300; -300];
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);


% Legenda colori
Legend_mag = ["A_d"; "A_n"; "\omega_{c,min}"; "G(j\omega)"];
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(GG_e,{omega_plot_min,omega_plot_max});
grid on; zoom on;


% Specifiche sovraelongazione (margine di fase)
omega_c_min = omega_Ta_MAX;
omega_c_MAX = omega_n_min;

phi_spec = Mf_spec - 180;
phi_low = -270; % lower bound per il plot

Bnd_Mf_x = [omega_c_min; omega_c_MAX; omega_c_MAX; omega_c_min];
Bnd_Mf_y = [phi_spec; phi_spec; phi_low; phi_low];
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;

% Legenda colori
Legend_arg = ["G(j\omega)"; "M_f"];
legend(Legend_arg);

% STOP qui per le specifiche



%% Rete anticipatrice

Mf_star = Mf_spec+1; % Mf_star = 69+1
%omega_c_star = omega_Ta_MAX; %circa 7,25
omega_c_star = 10;
[mar_omega_c_star, arg_omega_c_star, ~] = bode(GG_e, omega_c_star);

display(mar_omega_c_star);

mar_omega_c_star_db = 20*log10(mar_omega_c_star);

M_star = 10^(-mar_omega_c_star_db/20);


display(M_star);
phi_star = Mf_star - 180 - arg_omega_c_star;

%% Formule di inversione (rete anticipatrice)
% pi/180 per convertire da gradi a radianti o cosd/sind per usare i gradi
alpha_tau = (cosd(phi_star) - 1/M_star)/(omega_c_star*sind(phi_star));
tau = (M_star - cosd(phi_star))/(omega_c_star*sind(phi_star));
display(tau);
alpha = alpha_tau / tau ;
display(alpha);


% per avere alpha>0 si deve avere cos(phi_star)>1/M_star
check_flag = cosd(phi_star) - 1/(M_star);
display(cosd(phi_star));
display(1/M_star);
if check_flag < 0
    disp('Errore: alpha negativo');
    return;
end



%% Diagrammi di Bode con specifiche includendo regolatore dinamico

R_d = (1 + tau*s)/((1 + alpha*tau*s)); % rete anticipatrice (scenario di tipo B)

display(R_d);

LL = R_d*GG_e; % funzione di anello (incluso polo extra alte frequenze in GG_e)

figure(4);
%bode(LL,{omega_plot_min,omega_plot_max});
hold on;


% Specifiche su ampiezza
patch(Bnd_d_x, Bnd_d_y,'r','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_n_x, Bnd_n_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
patch(Bnd_Ta_x, Bnd_Ta_y,'b','FaceAlpha',0.2,'EdgeAlpha',0);
legend(Legend_mag);

% Plot Bode con margini di stabilità
margin(LL,{omega_plot_min,omega_plot_max});
grid on; zoom on;

% Specifiche su fase
patch(Bnd_Mf_x, Bnd_Mf_y,'g','FaceAlpha',0.2,'EdgeAlpha',0);
hold on;
legend(Legend_arg);

% STOP qui per sistema con controllore dinamico + specifiche




%% Check prestazioni in anello chiuso

% Funzione di sensitività complementare
FF = LL/(1+LL);

% gradino
WW = -25;

% Risposta al gradino
figure(5);

T_simulation = 5;
[y_step,t_step] = step(WW*FF, 5);
plot(t_step,y_step,'b');
grid on, zoom on, hold on;

% vincolo sovraelongazione
patch([0,T_simulation,T_simulation,0],[WW*(1+S_100_spec),WW*(1+S_100_spec),WW-5,WW-5],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
ylim([WW-5,0]);

% vincolo tempo di assestamento all'5%
LV = (evalfr(WW*FF,0)); % valore limite gradino: W*F(0)

patch([T_a5_spec,T_simulation,T_simulation,T_a5_spec],[WW*(1+0.05),WW*(1+0.05),WW-1,WW-1],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a5_spec,T_simulation,T_simulation,T_a5_spec],[WW*(1-0.05),WW*(1-0.05),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

Legend_step = ["Risposta al gradino"; "Vincolo sovraelongazione"; "Vincolo tempo di assestamento"];
legend(Legend_step);

% Diagramma di bode sistema in anello chiuso
figure(6);
bode(FF)
grid on, zoom on


%% Simulink point 4 - 5 and additional point 2 - 3

open("prova_simulink");
