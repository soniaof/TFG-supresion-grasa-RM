% test_fat_sat.m
% añadir directorios

codedir = cd;   % directorio donde está el código

datdir = '/Users/sonita/Desktop/4º/TFG/Datos/';

dirbloch = '/Users/sonita/Desktop/4º/TFG/Código TFG/Código CHESS/BlochSim/';
addpath(genpath(dirbloch));

dirirt = '/Users/sonita/Desktop/4º/TFG/Código TFG/Código CHESS/fessler/irt/';
cd(dirirt); 
setup; 
cd(codedir);
%% Crear parámetros para simulacion 

% variables que pertenecen at 1.5T
Sys=sys_sola();                             % carga parámetros del sistema, Siemens Sola Fit 1.5T
T1=t1_values_1_5;                           % carga valores del T1 a 1.5T
T2=t2_values_1_5;                           % carga valores del T2 a 1.5T
% pulso RF
Trf=5e-3;                                   % tiempo del pulso RF [sec], el pulso dura 5ms
time=[0:Sys.dt:Trf];                        % puntos del tiempo de pulse RF [sec],
% crea todos los instantes del pulso: 0, dt, 2dt… hasta 5 ms.

Npts=length(time);                          % longitud de pulso
flip=90;                                    % angulo de giro [deg]
b1amp=(flip*pi/180)/(2*pi*Sys.gamma*Trf);   % amplitud de pulso de RF [***calcular unidades], ecuación de Larmor
% unidades de amplitud de B1: [T], ya que Sys.gamma está en Hz/T, Hz = 1/s,
% y Trf en s

b1rf=[0; b1amp*ones((Npts-2),1); 0];          % crear pulso RF "hard" (rectangular), empieza y acaba en 0
gx=zeros(Npts,1); gy=gx; gz=gx;               % sin gradientes para pulso "hard"

%% Representar pulsos y gradientes 
%primera gráfica: pulso b1rf, que es rectangular, representado como linea recta (constante) 
% hasta 5ms, donde la línea azul es la intensidad de B1 y la naranja la fase de B1
% Los gradientes son campos magnéticos que cambian con la posición, la
% segunda gráfica los representa, pero como aún no hay, es una línea recta

figure; 
subplot(2,1,1);
yyaxis left
plot(time,abs(b1rf)); axis([time(1) time(end) 0 round(max(abs(b1rf))*1e7)/1e7]);
ylabel('|B1| ([T])') 
yyaxis right
plot(time,angle(b1rf)); axis([time(1) time(end) -pi pi]);
xlabel('Time (sec)'); ylabel('<B1 (rad)');
subplot(2,1,2);
plot(time,gx); hold on; plot(time,gy); plot(time,gz);
hold off; xlabel('Time (sec)'); ylabel('G_{x,y,z} (mT/m)');

%% Simulador de Bloch
FOVx=25;                    % [cm] 25 cm de espacio, distancia entre puntos: 25/128 = 0,195 cm
Nx=128;
x=[-Nx/2:Nx/2-1]*FOVx/Nx; x=x(:); %Esto crea muchos puntos a lo largo de un eje (de -12.5 cm a +12.5 cm aprox),
% es como si simulo muchos mini-protones colocados en distintas posiciones.
freq=zeros(size(x));        % [Hz], todos esos protones tienen frecuencia 0 Hz (esto es “on-resonance” con agua).
sens=ones([1 size(x)]);     % senisbilidad de RF, no lo utilizamos

% ***ojo: para correr el simulador de Bloch, hay que convertir las unidades 1T = 1e4 G ***

b1rfsim = b1rf * 1e4; %se convierte en gauss al multiplicar por 1e4
gxsim = gx * 0.1;
gysim = gy * 0.1;
gzsim = gz * 0.1; % (mT/m) -> (G/cm)  porque 1 mT/m = 0.1 G/cm

T1sim = T1.fat / 1000;
T2sim = T2.fat / 1000; %(elegir un tejido, por ejemplo grasa), divido entre 1000 para que sea en s
[Mx,My,Mz] = blochCim(b1rfsim,[gxsim, gysim, gzsim],Sys.dt,T1sim,T2sim,freq,x,0,sens);
Mxy=Mx+1i*My; %magnetización transversal (lo que “da señal”)

%% Representar simulación de Bloch "on-resonance" (a la frecuencia de agua), 0 Hz
% la primera gráfica es Mxy, donde el eje x es la posición espacial x(cm)
% de -12,5 a 12,5 cm aprox y el eje y (color azul) es la magnitud del eje transversal ,
% que es aprox 1,5 en la gráfica, esto significa que el pulso de RF ha creado magnetización transversal
% el color naranja es la fase.
%la segunda gráfica es Mz, la magnitud longitudinal que es 0, antes del pulso de RF era 1,
% en Mxy al revés, antes era 0 y ahora no.

figure;
subplot(2,1,1);
yyaxis left;
plot(x,abs(Mxy)); axis([x(1) x(end) -1 1]);
ylabel('|M_{xy}|'); title('M_{xy}');
yyaxis right;
plot(x,angle(Mxy)); axis([x(1) x(end) -pi pi]);
xlabel('X (cm)'); ylabel('<M_{xy}');
subplot(2,1,2);
plot(x,Mz); axis([x(1) x(end) -1 1]);
xlabel('X (cm)'); ylabel('|M_{z}|'); title('M_{z}');

%% **repetir la simulación del pulso hard con la frecuencia de grasa (calcularla desde el ppm en sys_sola.m)**

freqf = -Sys.gamma * Sys.field * Sys.fat * 1e-6; %sacamos del sys_sola el dato de Sys.fat=3.5 ppm y lo
% multiplicamos por sys.gamma=4.2576*1e7 y por Sys.field=1.5 al estar en
% campo de 1,5 T. Da 223,5240 Hz

freq_fat = freqf * ones(size(x));    % hace un vector con esa misma frecuencia para todos los protones.
[Mxf,Myf,Mzf] = blochCim(b1rfsim,[gxsim, gysim, gzsim],Sys.dt,T1sim,T2sim,freq_fat,x,0,sens);
Mxyf = Mxf + 1i*Myf;

figure;
subplot(2,1,1);
yyaxis left;
plot(x,abs(Mxyf)); axis([x(1) x(end) -1 1]);
ylabel('|M_{xy}|'); title('M_{xy}');
yyaxis right;
plot(x,angle(Mxyf)); axis([x(1) x(end) -pi pi]);
xlabel('X (cm)'); ylabel('<M_{xy}');
subplot(2,1,2);
plot(x,Mzf); axis([x(1) x(end) -1 1]);
xlabel('X (cm)'); ylabel('|M_{z}|'); title('M_{z}');
%% **cambiar el pulso para saturar la frecuencia de grasa y simular **

% frecuencia de la grasa (Hz) 
freqf = -Sys.gamma * Sys.field * Sys.fat * 1e-6;   % ~ -223.5 Hz a 1.5T

% se cambia la frecuencia del pulso RF para que esté centrado en la grasa(cambia la portadora),
% esto se hace multiplicando el pulso por una exponencial
b1rfsim_fat = b1rfsim .* exp(-1i*2*pi*freqf*time(:));

% simular grasa (su offset sigue siendo freqf)
% x contiene todas las posiciones espaciales donde está simulando
% protones, se crea un vector del mismo tamaño que x
%así todos los protones simulados tienen la misma frecuencia (223,5 Hz aprox)
freq_fat = freqf * ones(size(x)); 

[Mxf,Myf,Mzf] = blochCim(b1rfsim_fat,[gxsim,gysim,gzsim],Sys.dt,T1sim,T2sim,freq_fat,x,0,sens);
%si hacemos plot sale Mz en 0, ya que la grasa se ha saturado
% Mxyf = Mxf + 1i*Myf;
% 
% figure;
% subplot(2,1,1);
% yyaxis left;
% plot(x,abs(Mxyf)); 
% axis([x(1) x(end) 0 1]);
% ylabel('|M_{xy}|');
% title('M_{xy} de la grasa tras el pulso CHESS');
% 
% yyaxis right;
% plot(x,angle(Mxyf)); 
% axis([x(1) x(end) -pi pi]);
% xlabel('X (cm)');
% ylabel('fase M_{xy}');
% 
% subplot(2,1,2);
% plot(x,Mzf); 
% axis([x(1) x(end) -1 1]);
% xlabel('X (cm)');
% ylabel('M_z');
% title('M_z de la grasa tras el pulso CHESS');
% 
% idx0 = floor(Nx/2) + 1;
% fprintf('\n=== Grasa con pulso CHESS ===\n');
% fprintf('Frecuencia grasa = %.2f Hz\n', freqf);
% fprintf('|Mxy| en el centro = %.4f\n', abs(Mxyf(idx0)));
% fprintf('Mz en el centro = %.4f\n', Mzf(idx0));
%% **crear una simulación multi-epsectral para ver la saturación de frecuencias**

frange=[-300:10:100];                                   % rango de frecuencias [Hz] que vamos a probar de 10 en 10
% cubriendo las frecuencias de la grasa y el agua
Nf=length(frange);                                      % número de frecuencias
Mxspec=zeros(Nx,Nf); Myspec=Mxspec; Mzspec=Mxspec;      % Marcador de posición
Mxyspec=Mxspec;

for ii=1:Nf
    freq = frange(ii) * ones(size(x)); %ejemplo, si frange(ii) = -220 entonces freq = [-220 -220 -220 ...]

    % blochCim simula cómo evoluciona la magnetización de los espines al aplicar el pulso RF
    % para una determinada frecuencia (frange(ii)). Se repite la simulación para distintas frecuencias
    % para obtener la respuesta espectral del pulso y ver qué frecuencias quedan saturadas.
    [Mxspec(:,ii),Myspec(:,ii),Mzspec(:,ii)] = blochCim(... 
        b1rfsim_fat,[gxsim,gysim,gzsim],Sys.dt,T1sim,T2sim,freq,x,0,sens); 

    Mxyspec(:,ii)=Mxspec(:,ii)+1i*Myspec(:,ii); % construye Mxy = Mx + iMy, que es la señal detectable en MRI.
end

%con el for se muestra cómo responde el pulso según la frecuencia del spin.
Mxyspec_plot=squeeze(Mxyspec(end/2+1,:)); Mzspec_plot=squeeze(Mzspec(end/2+1,:));
figure;
subplot(2,1,1);
yyaxis left;
plot(frange,abs(Mxyspec_plot)); axis([frange(1) frange(end) -1 1]);
ylabel('|M_{xyf}|'); title('M_{xyf}');
yyaxis right;
plot(frange,angle(Mxyspec_plot)); axis([frange(1) frange(end) -pi pi]);
xlabel('freq (Hz)'); ylabel('<M_{xyf}');
subplot(2,1,2);
plot(frange,Mzspec_plot); axis([frange(1) frange(end) -1 1]);
xlabel('freq (Hz)'); ylabel('|M_{z}|'); title('M_{z}');

%% **crear una simulación temporal para ver la evolución de magnetización con la aplicación de RF**

Nt = length(time);                 % número de muestras temporales del pulso, time tiene los instantes desde 0 hasta Trf (5 ms)

%Creamos 3 vectores vacíos (llenos de ceros) para guardar, en cada tiempo: 
% Mx_t(kk) = Mx al tiempo kk, My_t (KK) = My al tiempo kk...
Mx_t = zeros(1,Nt); 
My_t = zeros(1, Nt);
Mz_t = zeros(1, Nt);

% (elegimos un solo punto espacial en vez de 128, por ejemplo el centro) se simula 1 solo spin en la posición x=0 cm
x0 = 0;
sens0 = 1;

% frecuencia (por ejemplo grasa centrada -> on-resonance)
freq0 = freqf;   

for kk = 1:Nt %Repite desde el primer punto de tiempo hasta el último.
    %Si kk=1: b1_tmp tiene solo el primer punto del pulso.
    %Si kk=10: b1_tmp tiene los 10 primeros puntos del pulso.
    %Si kk=Nt: b1_tmp tiene todo el pulso completo.
    b1_tmp = b1rfsim_fat(1:kk); %nos quedamos con el trozo de pulso desde 1 hasta kk
    g_tmp  = [gxsim(1:kk), gysim(1:kk), gzsim(1:kk)]; %igual con gradientes aunque da igual porque son 0

    [Mx1,My1,Mz1] = blochCim(b1_tmp, g_tmp, Sys.dt, T1sim, T2sim, freq0, x0, 0, sens0); %

    Mx_t(kk) = Mx1; %Guarda lo obtenido en el instante kk
    My_t(kk) = My1;
    Mz_t(kk) = Mz1;
end

Mxy_t = Mx_t + 1i*My_t;

figure;
plot(time, Mz_t, 'LineWidth', 2); grid on;
xlabel('Tiempo (s)'); ylabel('M_z');
title('Evolución temporal de M_z durante el pulso RF');
ylim([-1 1]);

figure;
plot(time, abs(Mxy_t), 'LineWidth', 2); grid on;
xlabel('Tiempo (s)'); ylabel('|M_{xy}|');
title('Evolución temporal de |M_{xy}| durante el pulso RF');
ylim([0 1]);

% al principio |Mxy| = 0 (no hay señal transversal)si el pulso actúa bien y |Mz| = 1
% luego |Mxy| sube, es decir crea señal y |Mz| baja


%% 
