% test_2D_spec.m
% añadir directorios
codedir=cd;                 % directorio donde está el código
datdir='/Users/sonita/Desktop/4º/TFG/Datos/';
dirbloch='/Users/sonita/Desktop/4º/TFG/Código TFG/Código SPSP/BlochSim/'; % directorio donde está el simulador de Bloch
addpath(genpath(dirbloch));
dirirt='/Users/sonita/Desktop/4º/TFG/Código TFG/Código SPSP/fessler/irt/';   % directorio donde está el toolbox IRT
cd(dirirt); setup; cd(codedir);
%% valores generales
Sys=sys_sola();                             % carga parámetros del sistema, Siemens Sola Fit 1.5T
T1=t1_values_1_5;                           % carga valores del T1 a 1.5T
T2=t2_values_1_5;                           % carga valores del T2 a 1.5T
FOVx=25;                                    % FOV [cm]
Nx=128;                                     % número de valores dimensión X
dx=FOVx/Nx;                                 % distancia entre vóxeles [cm]
x=[-Nx/2:Nx/2-1]*dx; x=x(:);                % vector de valores X [cm]
flip=30;                                    % angulo de giro [deg]
df=10;                                      % distancia entre frecuencias [Hz]
frange=[-300:df:100];                       % rango de frecuencias [Hz]
Nf=length(frange);                          % número de frecuencias

%% espacio k de excitacion (espacial, espectral)
Ninter=1;                                  % número de spiral interleaves
Fcoeff=4;                                  % coeficientes FOV
rmax=2;                                  % radius maximo de k-space [1/cm]

%Ninter (Número de interleaves): Se refiere a cuántas "vueltas" o "interleaves"
% realiza la espiral dentro del espacio k. Más interleaves significan una trayectoria
% más compleja y densa, lo que aumenta la selectividad espacial y mejora la precisión de la excitación en ciertas frecuencias.

%Fcoeff (Coeficiente de FOV): Controla el tamaño del campo de visión en el espacio espectral.
% Aumentar este valor aumenta el rango de frecuencias que se pueden explorar, lo cual
% es importante cuando se trata de separar la señal de diferentes tejidos (agua vs grasa).

%rmax (Radio máximo de k-space): Define el tamaño máximo de la trayectoria en el espacio k.
% Es una limitación que indica hasta dónde puede extenderse la trayectoria del espiral
% en términos de frecuencia o posición.

% Añadir las entradas a la función vds.m para los gradientes del sistema Sola
% [k_sp,g_sp,s_sp,time,r_sp,theta_sp] = vds(Sys.Smax*100,Sys.Gmax*0.1,Sys.dt,Ninter,Fcoeff,rmax);
% %smax= amplitud máxima de los gradiente, gmax= máxima fuerza de los gradientes.
% Nt=length(k_sp);
% kx=real(k_sp(end:-1:1));                    % espacio k espacial de dimensión X [1/cm] 
% ky=imag(k_sp(end:-1:1));                    % espacio k de dimensión Y, no lo utilizamos aquí [1/cm]
%kt=time(end:-1:1);                          % espacio k espectral [s]  % como se define espacio k espectral? ** mirar tesis de Sydney**

[kx1,ky1] = espiralArquimedes(0, 0.1, 5, 1000);
[kx2,ky2,theta,r] = espiralDensidadVariable(5,5,1000,2)
kx=kx2(end:-1:1); ky=ky2(end:-1:1);
g=k2g([kx; ky], FOVx);
gx=g(1,:); gy=g(2,:);
Nt=length(kx);
time=linspace(0,(Nt-1)*Sys.dt,Nt);
kt=time(end:-1:1);                          % espacio k espectral [s]  % como se define espacio k espectral? ** mirar tesis de Sydney**


%La unidad 1/cm para kx (el espacio k espacial en la dirección x) se debe a que el
% espacio k está relacionado con la frecuencia espacial, que está inversamente relacionada
% con la distancia.


%% establecer target, magnetización ideal

dtarget=zeros(Nf,Nx);                       % provisional 
dfat=[1:10, 35:75, 99:118];                 % valores espacial de grasa
dwater=[11:34,76:98,119:128];               % valores espacial de agua
freqband=15;                                % ancho de banda de inclusión de grasa/agua [Hz],
%se toman todas las frecuencias dentro de un rango de ±15 Hz de sus frecuencias centrales
zoneband=25;                                % ancho de banda de transición de no grasa/no agua [Hz]
%rango de frecuencias intermedias donde la señal está en transición,
% no está ni en 0Hz (agua) ni en -220 Hz (grasa)
%Estamos permitiendo un margen de ±25 Hz alrededor de las frecuencias de agua (0 Hz)
% y grasa (~-220 Hz), para hacer una transición suave.

fatfreq=-(Sys.fat*1e-6*Sys.gamma*Sys.field);     % frecuencia de grasa calculado para 3T, definirla negativa al respecto de agua [Hz] 
fatval=find(abs(frange-fatfreq)<freqband);  % índices de frecuencias de grasa
watval=find(abs(frange-0)<freqband);        % índices de frecuencias de agua
roival=find((abs(frange-0)<zoneband) | (abs(frange-fatfreq)<zoneband)); % índices de inclusión en el ROI para el diseño
dtarget(watval,dwater)=0;                   % angulo de giro deseado para agua [rad]
dtarget(fatval,dfat)=flip*pi/180;                % angulo de giro desado para grasa **saturación** [rad]
roi=logical(zeros(size(dtarget)));          % definir ROI
roi(roival,:)=true;

figure;
subplot(1,2,1);
imagesc(x,frange,dtarget); xlabel('X [cm]'); ylabel('frec. [Hz]');
h=colorbar; ylabel(h,'[rad]'); title('Target, Saturación de Grasa');
subplot(1,2,2);
imagesc(x,frange,roi); colormap(gca,'gray'); xlabel('X [cm]'); ylabel('frec. [Hz]');
title('ROI, Frecuencias de Inclusión');
% ** explicar lo que se está mostrando en estas figuras **

% La figura izquierda muestra el target de magnetización en el dominio
% espectral-espacial (X vs frecuencia). Se observa que solo se impone un
% ángulo de giro distinto de cero (~90°, 90º=2π/2​≈1.57 rad) en la banda de la grasa (~−220 Hz)
% y en ciertas posiciones espaciales, mientras que el agua (0 Hz) no se excita.
%
% La figura derecha muestra la ROI (región de interés), donde se incluyen
% únicamente las bandas de frecuencia del agua (0 Hz) y de la grasa (~−220 Hz).
% El diseño del pulso se optimiza solo en estas regiones.
%% Crear matriz de sistema, A, optimizar pulso espectral-espacial
% para esta parte, tendrías que revisar las secciones 2.2 y 2.3 de la tesis
% de Sydney

m0=1;                                       % definimos la magnetización empezando longitudinal
[F,X]=ndgrid(frange,x);                     % replicamos valores par n>1 dimensiones del diseño
r_spsp=[F(:) X(:)];                         % dimensiones r de la matriz del sistema, A
k_spsp=[kt; kx];                            % dimensiones k de la matriz del sistema, A
A = (1i*2*pi*4258*m0*Sys.dt)*exp(1i*2*pi*(r_spsp*k_spsp));   % añadir lo que falta en la matriz, la parte que no está en la exponencial
% % El factor fuera de la exponencial proviene de la aproximación STA:
% i * gamma * m0 * dt
%
% donde:
% - 4258 corresponde a gamma (constante giromagnética) en Hz/G -> 42576000 Hz/T
% - Sys.dt es el paso temporal (Δt)
% - 1i es el factor complejo de la aproximación STA (small tip angle)
% - 2*pi aparece porque se trabaja en frecuencia en Hz (no en rad/s)

% antes de seguir, ¿qué tamaño y cuáles unidades tiene A? ¿tiene sentido?
% A tiene tamaño (Nf*Nx) x Nt. Esto es porque r_spsp tiene tamaño
% (Nf*Nx) x 2 y k_spsp tiene tamaño 2 x Nt, luego r_spsp*k_spsp es
% (Nf*Nx) x Nt.
%
% En cuanto a unidades: F está en [Hz], kt en [s], X en [cm] y kx en [1/cm],
% de modo que F*kt y X*kx son adimensionales y la exponencial tiene
% sentido. La A es 1/G, el pulso b está en G.
% El factor fuera de la exponencial corresponde a la aproximación STA.

beta=1e-1;              % parámetro de regularizador, formula 2.13 tesis, la W es el ROI
% solución "closed-form" para weighted least-squares
bopt=pinv(A'*diag(roi(:))*A+beta*eye(Nt))*A'*diag(roi(:))*dtarget(:); % closed for solution for regularized WLS [V/rad]*[rad/V]*[V/rad]*rad = V
% aproximación STA de magnetización
mopt=reshape(A*bopt,[Nf Nx]);
% función anónima para calcular NRMSE
func.nrmse=@(m,d,r) norm(abs(d(r))-abs(m(r)))/norm(abs(d(r)));
nrmse_opt=func.nrmse(mopt,dtarget,roi);


%% Figuras
figure;
% ** crea una figura mostrando la magnitud y el fase de pulso y los espacios
% k espectral y esapcial por tiempo. Incluir unidades **

subplot(2,2,1);
plot(time*1e3,abs(bopt),'LineWidth',1.5);
xlabel('Tiempo [ms]');
ylabel('|b(t)|');
title('Magnitud del pulso RF');

subplot(2,2,2);
plot(time*1e3,angle(bopt),'LineWidth',1.5);
xlabel('Tiempo [ms]');
ylabel('Fase [rad]');
title('Fase del pulso RF');

subplot(2,2,3);
plot(time*1e3,kx,'LineWidth',1.5);
xlabel('Tiempo [ms]');
ylabel('k_x [1/cm]');
title('k-space espacial');

subplot(2,2,4);
plot(time*1e3,kt,'LineWidth',1.5);
xlabel('Tiempo [ms]');
ylabel('k_t [s]');
title('k-space espectral');

% La figura muestra el pulso RF obtenido. La magnitud indica la intensidad
% del pulso en el tiempo y la fase muestra cómo cambia su comportamiento
% complejo.
%
% Las gráficas de kx y kt representan cómo se recorre el espacio (kx) y la
% frecuencia (kt) a lo largo del tiempo. kx varía de forma suave debido a
% las limitaciones de los gradientes, mientras que kt crece de forma lineal
% porque depende directamente del tiempo.
%
% En las imágenes, el color representa el ángulo de giro: azul indica poca
% excitación y amarillo excitación alta. Se observa que el resultado no
% coincide perfectamente con el target, especialmente en la región de
% interés, lo que indica que el diseño del pulso no es ideal.

figure; 
subplot(1,2,1);
imagesc(x,frange,dtarget,[0 (flip*pi/180)]); xlabel('X [cm]'); ylabel('frec. [Hz]');
h=colorbar; ylabel(h,'[rad]'); title('Target, Saturación de Grasa');
% subplot(1,3,2)
% imagesc(x,frange,abs(mopt),[0 (flip*pi/180)]); xlabel('X [cm]'); ylabel('frec. [Hz]');
% h=colorbar; ylabel(h,'[rad]'); title('Magnetización, Toda');
subplot(1,2,2);
imagesc(x,frange,abs(mopt).*roi,[0 (flip*pi/180)]); xlabel('X [cm]'); ylabel('frec. [Hz]');
h=colorbar; ylabel(h,'[rad]'); title(sprintf('Magnetización, con ROI, NRMSE: %4.3f',nrmse_opt));

%% Próximo pasos

%target → lo que quieres
%mopt → lo que sale
%ROI → dónde miras
% El target define la magnetización deseada en función de la posición y la
% frecuencia. En este caso, se quiere excitar la grasa (ángulo de giro ~90°)
% y no excitar el agua.

% El código construye un modelo del sistema (matriz A) y calcula el pulso RF
% óptimo (bopt) que mejor reproduce el target. Después se obtiene la
% magnetización resultante (mopt) y se compara con el target mediante el
% NRMSE.
        
% 1. iterar por varios valores de beta y crea figuras y observaciones en el
% valor de beta y el efecto en el NRMSE y la magnitude del pulso

% beta controla el compromiso entre ajuste al target y suavidad del pulso.
% Valores pequeños ajustan mejor pero pueden dar pulsos más intensos,
% mientras que valores grandes suavizan el pulso pero empeoran el ajuste.

%bopt=pinv(A'*diag(roi(:))*A+beta*eye(Nt))*A'*diag(roi(:))*dtarget(:);

% cuanto menor es beta cada vez menos magnitud del pulso y NRMSE y cada vez mas
% norma de b (norm(bopt))


betas = [1e-4, 1e-2, 1e-1, 1];

figure('Position', [100 100 1400 900]);
sgtitle('Efecto de \beta en el diseño del pulso SPSP');

for i = 1:length(betas)
    b = betas(i);
    
    bopt_i = pinv(A'*diag(roi(:))*A + b*eye(Nt)) * A'*diag(roi(:))*dtarget(:);
    mopt_i = reshape(A*bopt_i, [Nf Nx]);
    nrmse_i = func.nrmse(mopt_i, dtarget, roi);

    fprintf('beta=%.0e | max(abs(bopt))=%.4f | norm(bopt)=%.4f | NRMSE=%.4f\n', b, max(abs(bopt_i)), norm(bopt_i), nrmse_i);

    subplot(4, 3, (i-1)*3 + 1);
    imagesc(x, frange, dtarget, [0 (flip*pi/180)]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]');
    h=colorbar; ylabel(h,'[rad]');
    title(sprintf('Target | \\beta=%.0e', b));

    subplot(4, 3, (i-1)*3 + 2);
    imagesc(x, frange, abs(mopt_i), [0 (flip*pi/180)]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]');
    h=colorbar; ylabel(h,'[rad]');
    title(sprintf('Magnetización | \\beta=%.0e', b));

    subplot(4, 3, (i-1)*3 + 3);
    imagesc(x, frange, abs(mopt_i).*roi, [0 (flip*pi/180)]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]');
    h=colorbar; ylabel(h,'[rad]');
    title(sprintf('ROI | NRMSE=%.3f | ||b||=%.3f', nrmse_i, norm(bopt_i)));
end

% En este apartado se estudia el efecto del parámetro de regularización beta
% sobre el diseño del pulso espectral-espacial.
%
% Se prueban 4 valores: beta = 1e-4, 1e-2, 1e-1, 1
%
% Para cada valor se recalcula el pulso óptimo (bopt) y la magnetización
% resultante usando la aproximación lineal m=A*b (small tip angle, STA).
% Se muestran el target, la magnetización completa y la magnetización
% dentro de la ROI.
%
% Se calculan dos métricas por cada beta:
% - NRMSE: mide qué tan bien se ajusta la magnetización al target dentro
%   de la ROI. Cuanto menor, mejor ajuste.
% - ||b||: norma del pulso, mide la energía total del pulso RF.
%   Cuanto menor, más suave es el pulso.
%
% Resultados esperados:
% - Beta pequeño: NRMSE bajo (buen ajuste) pero ||b|| alto (pulso intenso).
% - Beta grande: NRMSE alto (peor ajuste) pero ||b|| bajo (pulso suave).
% - Beta controla el compromiso entre calidad del diseño y energía del pulso.

%% 

% 2. cambia el dtarget/ROI para incluir mas o menos frecuencias y con menos
% giro de angulo y observar los resultados con la optimización

%dtarget(fatval,dfat)=flip*pi/180;  

%roival=find((abs(frange-0)<zoneband) | (abs(frange-fatfreq)<zoneband)); % índices de inclusión en el ROI para el diseño


% Al reducir el parámetro zoneband, la ROI se hace más estrecha en
% frecuencia, por lo que el problema de optimización es más fácil. En este
% caso se observa una mayor intensidad de la magnetización (zonas más
% amarillas), ya que el pulso puede concentrarse mejor en las regiones de
% interés.

% Al aumentar zoneband, la ROI se amplía y el pulso debe cumplir el target
% en un rango mayor de frecuencias, lo que dificulta el diseño y reduce la
% calidad de la excitación.


%% 2. Cambiar dtarget/ROI y simular con Bloch

% gradientes del espiral calculados a partir de la trayectoria actual
gxsim = gx(:);
gysim = gy(:);
gzsim = zeros(size(gxsim));

% tejido: grasa
T1sim = T1.fat / 1000;
T2sim = T2.fat / 1000;
sens  = ones([1 length(x)]);

% casos a probar: [flip]
flips = [30, 60, 90, 120];  % Solo flip cambia
zoneband = 25;  % zoneband no cambia

figure('Position',[100 100 1400 900]);
sgtitle('Efecto de flip - Simulación Bloch');

for i = 1:length(flips)

    % recalcular target y ROI
    dtarget_i = zeros(Nf, Nx);
    dtarget_i(fatval, dfat) = flips(i) * pi/180;

    % ROI fija (sin cambios en zoneband)
    roival_i = find(abs(frange - 0) < zoneband | abs(frange - fatfreq) < zoneband);
    roi_i = logical(zeros(size(dtarget_i)));
    roi_i(roival_i, :) = true;

    % recalcular pulso óptimo
    bopt_i = pinv(A'*diag(roi_i(:))*A + beta*eye(Nt)) * A'*diag(roi_i(:))*dtarget_i(:);

    % simular con Bloch para cada frecuencia
    Mxyspec = zeros(Nx, Nf);
    for ii = 1:Nf
        freq_ii = frange(ii) * ones(size(x));
        [Mx, My, ~] = blochCim(bopt_i, [gxsim, gysim, gzsim], Sys.dt, T1sim, T2sim, freq_ii, x, 0, sens);
        Mxyspec(:, ii) = Mx + 1i * My;
    end

    NRMSE_val = norm(abs(Mxyspec)'.*roi_i - dtarget_i) / norm(dtarget_i);
    fprintf('flip=%d | max(abs(bopt))=%.4f | norm(bopt)=%.4f | NRMSE=%.4f\n', flips(i), max(abs(bopt_i)), norm(bopt_i), NRMSE_val);

    % Target
    subplot(4, 3, (i - 1) * 3 + 1);
    imagesc(x, frange, dtarget_i, [0 flips(i) * pi / 180]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]'); colorbar;
    title(sprintf('Target | flip=%d | zoneband=%d', flips(i), zoneband));

    % Magnetización Bloch
    subplot(4, 3, (i - 1) * 3 + 2);
    imagesc(x, frange, abs(Mxyspec)', [0 flips(i) * pi / 180]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]'); colorbar;
    title('Magnetización Bloch');

    % ROI
    subplot(4, 3, (i - 1) * 3 + 3);
    imagesc(x, frange, abs(Mxyspec)'.*roi_i, [0 flips(i) * pi / 180]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]'); colorbar;
    NRMSE_val = norm(abs(Mxyspec)'.*roi_i - dtarget_i) / norm(dtarget_i);  % Calcular NRMSE
    title(sprintf('ROI | ||b||=%.3f | NRMSE=%.3f', norm(bopt_i), NRMSE_val));

end

% Este código analiza cómo cambia la respuesta de la técnica de supresión espectral
% cuando se varía el ángulo de giro (flip) de un pulso RF, mientras que el ancho
% de banda de la ROI (zoneband) se mantiene constante.
%
% Se prueba un conjunto de valores para el ángulo de giro (flip) en 30°, 60°, 
% 90° y 120° para observar cómo afecta la magnetización de la grasa.
% 
% Los pasos incluyen:
% 1. Recalcular el target y la ROI para cada valor de flip.
% 2. Calcular el pulso óptimo y la magnetización utilizando el simulador de Bloch.
% 3. Mostrar la magnetización resultante y calcular el NRMSE entre la magnetización 
%    obtenida y el target, reflejando el desempeño de la simulación.
% 
% Resultados esperados:
% - A medida que aumenta el ángulo de giro (flip), se observa un incremento en 
%   la norma del pulso óptimo (||b||), lo que indica un mayor esfuerzo en la excitación.
% - El NRMSE proporciona una métrica de la diferencia entre la magnetización deseada
%   y la resultante, lo que permite evaluar la precisión del pulso optimizado.
   
%% 

% 3. jugar con la función vds.m para crear un espiral más dinámico (más turnos) y
% observar resultados con la optimización (intentar no pasar una
% trayectoría más larga que 10 ms)

% Ninter es número de vueltas del espiral en k-space.
%A mayor Ninter, la trayectoria es más compleja y puede mejorar
% la selectividad, aunque aumenta la duración del pulso.

%vds(smax, gmax, T, N, Fcoeff, rmax)
%como smax, gmax y T son fijos puedo variar lo demás: N, Fcoeff, rmax
