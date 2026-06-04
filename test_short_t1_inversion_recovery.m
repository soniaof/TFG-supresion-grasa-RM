% test_short_t1_inversion_recovery.m
% añadir directorios
codedir=cd;                 % directorio donde está el código
datdir='/Users/sonita/Desktop/4º/TFG/Datos/';
dirbloch='/Users/sonita/Desktop/4º/TFG/Código TFG/Código STIR/BlochSim/'; % directorio donde está el simulador de Bloch
addpath(genpath(dirbloch));
dirirt='/Users/sonita/Desktop/4º/TFG/Código TFG/Código STIR/fessler/irt/';                      % directorio donde está el toolbox IRT
cd(dirirt); setup; cd(codedir);
%% Crear parámetros para simulacion 

% variables que pertenecen at 1.5T
Sys=sys_sola();                               % carga parámetros del sistema, Siemens Sola Fit 1.5T
T1=t1_values_1_5;                           % carga valores del T1 a 1.5T
T2=t2_values_1_5;                           % carga valores del T2 a 1.5T
% pulso RF
Trf=10e-3;                                  % tiempo del pulso RF [sec]
time=[0:Sys.dt:Trf];                        % puntos del tiempo de pulse RF [sec]
Npts=length(time);                          % longitud de pulso

%% Reproducir pulso adiabático que invertie la señal (pulso de ineversión)
% aquí, no vamos a simular el pulso después para STIR pero sirve para ver cómo funcionan los pulsos de inversión
b1max=0.2;                                  % magnitud del pulso de inversion [G]
b1bw=2e3;                                   % ancho de banda del pulso [Hz]
beta=1e3;                                   % parámetro 
doBloch=1;                                  % si/no hacer simulación
[rf_adi,freq_adi,phase_adi]=adiabatic(b1max,b1bw,beta,Trf,Sys.dt,1); % crea pulso de inversión adiabática (función por la Universidad de Stanford)

% aquí,¿cuál es la diferencia de fase y frecuencia con el pulso?
% La fase representa el ángulo instantáneo del campo RF en el plano transversal xy,
% mientras que la frecuencia instantánea describe la velocidad de cambio de esa fase en el tiempo.
% Matemáticamente, la frecuencia es la derivada temporal de la fase:
%
%           f(t) = (1/2π) dφ(t)/dt
%
% En los pulsos adiabáticos la frecuencia varía progresivamente (frequency sweep) alrededor de la resonancia,
% lo que produce un cambio continuo de fase y permite invertir la magnetización.

% Aunque la fase del pulso RF se define en el plano transversal (xy),
% ya que el campo de RF B1 actúa perpendicularmente al campo
% principal B0 (es decir en xy, plano transversal),
% este campo provoca la rotación del vector de magnetización completo.
% Como resultado, el pulso puede invertir la magnetización longitudinal (Mz → -Mz),
% que es la componente relevante en técnicas de inversión como STIR.

% FIGURA 1: Características del pulso adiabático Silver-Hoult

% Esta figura muestra las propiedades temporales del pulso RF adiabático generado.

% (1) Amplitud B1(t):
% Representa la envolvente del campo de radiofrecuencia aplicado. La amplitud
% aumenta suavemente hasta un máximo y luego disminuye, lo que ayuda a cumplir
% la condición adiabática y permite que la magnetización siga el campo efectivo.

% (2) Fase del pulso:
% Muestra la fase instantánea del campo RF. La representación aparece con
% discontinuidades debido al wrapping de fase entre -pi y pi, pero en realidad
% la fase evoluciona de forma continua durante el pulso.

% (3) Frecuencia instantánea:
% Corresponde a la derivada temporal de la fase. El pulso realiza un barrido
% de frecuencia desde valores negativos hasta positivos alrededor de la
% frecuencia de resonancia. Este barrido permite invertir la magnetización
% de manera robusta dentro de un rango de frecuencias.

% FIGURA 2: Respuesta de la magnetización simulada con las ecuaciones de Bloch

% Esta figura muestra el efecto del pulso adiabático sobre la magnetización
% para distintos offsets de frecuencia. Se representan las tres componentes
% de la magnetización final tras aplicar el pulso:

% Mx (azul)  -> componente transversal en el eje x
% My (verde) -> componente transversal en el eje y
% Mz (rojo)  -> componente longitudinal

% La componente longitudinal Mz muestra que la magnetización se invierte
% (Mz ≈ -1) dentro de un rango de frecuencias alrededor de la resonancia,
% mientras que fuera de ese rango permanece aproximadamente en +1.

% Esto demuestra que el pulso funciona como un pulso de inversión selectivo
% en frecuencia. Las oscilaciones en Mx y My corresponden a la rotación de
% la magnetización en el plano transversal durante el proceso de inversión.

%% Calculación de recuperación de T1

time_sim=[0:Sys.dt:3];                      % puntos del tiempo simulación de recuperación T1 [sec], es de 0 a 3 sec

T1 = t1_values_1_5;                      % cargar valores T1 (ms)

% convertir a segundos
T1_fat = T1.fat/1000;
T1_muscle = T1.muscle/1000;
T1_water = T1.csf/1000;

Mz_fat = 1 - 2*exp(-time_sim/T1_fat);                              % buscar formula para recuperación T1 utilizando el vector de tiempo time_sim y valor de T1 de grasa *cuidado con unidades*
Mz_muscle = 1 - 2*exp(-time_sim/T1_muscle);                        % repetir con músculo
Mz_water = 1 - 2*exp(-time_sim/T1_water);                          % repetir con agua (utiliza valor de CSF)


% Se utilizan valores de T1 obtenidos de la función t1_values_1_5 (en ms).
% Estos se convierten a segundos para la simulación. La recuperación longitudinal tras un pulso de inversión
% sigue la ecuación:

%       Mz(t) = 1 - 2exp(-t/T1)

% lo que permite simular la recuperación de diferentes tejidos como grasa,
% músculo y agua (CSF) en función de sus tiempos de relajación T1.
%% Representar pulsos y gradientes
figure; 
plot(time_sim,Mz_fat,'b'); hold on;
plot(time_sim,Mz_muscle,'Color',[0.25 1 0]);
plot(time_sim,Mz_water,'Color',[0.5 0.5 0]);
plot(time_sim,zeros(size(time_sim)),'r--');
legend('Grasa','Músculo','Agua');
xlabel('Tiempo [s]');
ylabel('Magnetización longitudinal M_z');

% ¿está correcto la linea de recuperación de T1? (i.e., empieza con la
% longitud de magnetización Mz en -1 y llega a 1 al final) ¿cómo puedes
% cambiar la formula para que sea correcto?

% No, inicialmente la línea de recuperación no es correcta porque los valores de T1 estaban en milisegundos,
% mientras que el vector de tiempo time_sim está en segundos. Para que la recuperación sea correcta, hay que
% convertir T1 a segundos y usar esos valores en la fórmula:
% Mz(t) = 1 - 2*exp(-t/T1 en segundos).
% Así, la magnetización empieza en Mz = -1 y tiende a Mz = +1 con el tiempo.

% cuándo esté corregido, con el ratón, en la figura, determina los valores de tiempo que
% anula grasa, músculo, y agua. Añade esta información a la figura en una
% caja de texto o en el título

% Los tiempos de anulación obtenidos de la figura son aproximadamente:
% grasa ≈ 0.18 s = 180 ms
% músculo ≈ 0.60 s
% agua ≈ 1.7 s

txt = sprintf('TI grasa \\approx 0.18 s\nTI músculo \\approx 0.60 s\nTI agua \\approx 1.7 s');

text(1.6,-0.35,txt, ...
    'BackgroundColor','white', ...
    'EdgeColor','black', ...
    'Interpreter','tex');

Fat0 = find(min(abs(Mz_fat))==abs(Mz_fat));
Water0 = find(min(abs(Mz_water))==abs(Mz_water));
Muscle0 = find(min(abs(Mz_muscle))==abs(Mz_muscle));