
% Extensión de test_2D_spec añadiendo dimensión espacial Y

%% Añadir directorios
codedir = cd;
dirbloch = '/Users/sonita/Desktop/4º/TFG/Código TFG/Código SPSP/BlochSim/';
addpath(genpath(dirbloch));
dirirt = '/Users/sonita/Desktop/4º/TFG/Código TFG/Código SPSP/fessler/irt/';
cd(dirirt); setup; cd(codedir);

%% Parámetros generales
Sys    = sys_sola();
T1     = t1_values_1_5;
T2     = t2_values_1_5;

FOVx = 25;  Nx = 32;  dx = FOVx/Nx;
x = [-Nx/2 : Nx/2-1] * dx;  x = x(:);   % vector X [cm]

% NUEVO: dimensión Y igual que X
FOVy = 25;  Ny = 32;  dy = FOVy/Ny;
y = [-Ny/2 : Ny/2-1] * dy;  y = y(:);   % vector Y [cm]

flip   = 30;
df     = 10;
frange = [-300:df:100];
Nf     = length(frange);
m0     = 1;

fprintf('Dimensiones: Nf=%d, Nx=%d, Ny=%d\n', Nf, Nx, Ny);
fprintf('Tamaño de A: (%d x Nt)\n', Nf*Nx*Ny);

%% Trayectoria espiral k-space
Ninter = 1;
Fcoeff = 4;
rmax   = 2;

[k_sp, g_sp, s_sp, time, r_sp, theta_sp] = vds(Sys.Smax*100, Sys.Gmax*0.1, ...
    Sys.dt, Ninter, Fcoeff, rmax);

Nt = length(k_sp);
kx = real(k_sp(end:-1:1));   % k-space espacial X [1/cm]
ky = imag(k_sp(end:-1:1));   % k-space espacial Y [1/cm] <- NUEVO
kt = time(end:-1:1);          % k-space espectral [s]
%Antes solo usaba kx, ahora uso también ky porque tengo la dimensión Y.
% kt es el tiempo invertido, que actúa como dimensión espectral.
% El end:-1:1 invierte el vector porque la espiral va de fuera hacia dentro y
% se necesita recorrerla al revés.
fprintf('Longitud del pulso: %d puntos = %.2f ms\n', Nt, time(end)*1000);

%% Establecer target y ROI
freqband = 15;
zoneband = 25;
fatfreq  = -(Sys.fat * 1e-6 * Sys.gamma * Sys.field);

fatval = find(abs(frange - fatfreq) < freqband);
watval = find(abs(frange - 0)       < freqband);
roival = find((abs(frange - fatfreq) < zoneband) | (abs(frange - 0) < zoneband));

% Zonas espaciales escaladas a Nx=32
dfat   = [1:3, 9:19, 25:30];
dwater = [4:8, 20:24, 31:32];

% dtarget 3D: (Nf x Nx x Ny)
% En 1D era: dtarget(fatval, dfat) = flip*pi/180
% En 2D es:  dtarget(fatval, dfat, :) = flip*pi/180
dtarget = zeros(Nf, Nx, Ny);
dtarget(fatval, dfat, :) = flip * pi/180; %en las frecuencias de grasa (fatval),
% en las columnas de grasa (dfat), y en TODAS las posiciones Y (:) pon el ángulo
% de flip. El resultado es que el target es uniforme en Y, es decir,
% las franjas de grasa son iguales para cualquier Y que mires.

% ROI 3D: (Nf x Nx x Ny)
% En 1D era: roi(roival, :) = true
% En 2D es:  roi(roival, :, :) = true
roi = false(Nf, Nx, Ny);
roi(roival, :, :) = true;

% Visualizar corte central en Y
iy_c = round(Ny/2);
figure('Name','Target y ROI - corte Y central');
subplot(1,2,1);
imagesc(x, frange, dtarget(:,:,iy_c), [0 flip*pi/180]);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title(sprintf('Target 2D, corte Y=%.1f cm', y(iy_c)));

subplot(1,2,2);
imagesc(x, frange, double(roi(:,:,iy_c)));
colormap(gca,'gray'); xlabel('X [cm]'); ylabel('Frec. [Hz]');
title('ROI 2D, corte Y central');

%% Construir matriz del sistema A
%
% CAMBIO CLAVE respecto al original:
%
%   ANTES (1D espacial):
%     [F, X]    = ndgrid(frange, x)
%     r_spsp    = [F(:) X(:)]          tamaño (Nf*Nx) x 2
%     k_spsp    = [kt; kx]             tamaño 2 x Nt
%
%   AHORA (2D espacial):
%     [F, X, Y] = ndgrid(frange, x, y)
%     r_spsp    = [F(:) X(:) Y(:)]     tamaño (Nf*Nx*Ny) x 3
%     k_spsp    = [kt; kx; ky]         tamaño 3 x Nt
%
% La exponencial sigue siendo adimensional:
%   F*kt [Hz*s] + X*kx [cm*1/cm] + Y*ky [cm*1/cm] -> sin unidades

fprintf('\nConstruyendo matriz A (%d x %d)...\n', Nf*Nx*Ny, Nt);

[F, X, Y] = ndgrid(frange, x, y);

r_spsp = [F(:), X(:), Y(:)];
k_spsp = [kt; kx; ky];

A = (1i * 2*pi * 4258 * m0 * Sys.dt) * exp(1i * 2*pi * (r_spsp * k_spsp));

fprintf('Matriz A construida: %d x %d\n', size(A,1), size(A,2));

%% Optimizar pulso
beta = 1e-1;

dtarget_vec = dtarget(:);
roi_vec     = roi(:);

fprintf('Optimizando pulso...\n');

% Sin crear matriz diag gigante (evita error de memoria)
Aw  = A .* roi_vec;
AtA = Aw' * A;
Atd = A'  * (roi_vec .* dtarget_vec);
bopt = (AtA + beta * eye(Nt)) \ Atd;

% Magnetización resultante (aproximación STA)
mopt = reshape(A * bopt, [Nf, Nx, Ny]); %calcula la magnetización resultante multiplicando
% la matriz A por el pulso, y luego la reorganiza de vector a array 3D

% NRMSE
func_nrmse = @(m,d,r) norm(abs(d(r)) - abs(m(r))) / norm(abs(d(r)));
nrmse_2D   = func_nrmse(mopt, dtarget, roi);

fprintf('NRMSE 2D: %.4f\n', nrmse_2D);
fprintf('Norma del pulso ||b||: %.4f\n', norm(bopt));

%% Figuras del pulso RF
figure('Name','Pulso RF optimizado');
subplot(2,2,1);
plot(time*1e3, abs(bopt), 'LineWidth', 1.5);
xlabel('Tiempo [ms]'); ylabel('|b(t)|');
title('Magnitud del pulso RF');

subplot(2,2,2);
plot(time*1e3, angle(bopt), 'LineWidth', 1.5);
xlabel('Tiempo [ms]'); ylabel('Fase [rad]');
title('Fase del pulso RF');

subplot(2,2,3);
plot(time*1e3, kx, 'b', 'LineWidth', 1.5); hold on;
plot(time*1e3, ky, 'r', 'LineWidth', 1.5); hold off;
xlabel('Tiempo [ms]'); ylabel('k [1/cm]');
legend('k_x','k_y'); title('k-space espacial X e Y');

subplot(2,2,4);
plot(time*1e3, kt, 'LineWidth', 1.5);
xlabel('Tiempo [ms]'); ylabel('k_t [s]');
title('k-space espectral');

%% Figura resultado - corte central en Y
figure('Name','Resultado SPSP 2D - corte Y central');
subplot(1,3,1);
imagesc(x, frange, dtarget(:,:,iy_c), [0 flip*pi/180]);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title('Target (corte Y central)');

subplot(1,3,2);
imagesc(x, frange, abs(mopt(:,:,iy_c)), [0 flip*pi/180]);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title('Magnetización (corte Y central)');

subplot(1,3,3);
imagesc(x, frange, abs(mopt(:,:,iy_c)) .* roi(:,:,iy_c), [0 flip*pi/180]);
xlabel('X [cm]'); ylabel('Frec. [Hz]');
h = colorbar; ylabel(h,'[rad]');
title(sprintf('ROI | NRMSE=%.3f | ||b||=%.3f', nrmse_2D, norm(bopt)));

%% Figura imagen espacial X-Y en banda de grasa
mopt_grasa    = squeeze(sum(abs(mopt(fatval,:,:)), 1));
dtarget_grasa = squeeze(sum(dtarget(fatval,:,:), 1));

figure('Name','Imagen espacial 2D - banda de grasa');
subplot(1,2,1);
imagesc(x, y, dtarget_grasa.'); colormap(gca,'parula'); axis image; colorbar;
xlabel('X [cm]'); ylabel('Y [cm]');
title('Target espacial X-Y (banda grasa)');

subplot(1,2,2);
imagesc(x, y, mopt_grasa.'); colormap(gca,'parula'); axis image; colorbar;
xlabel('X [cm]'); ylabel('Y [cm]');
title(sprintf('Magnetización X-Y (banda grasa) | NRMSE=%.3f', nrmse_2D));

%% VISUALIZACIÓN 3D - Vista 1: Varios cortes en Y
% Muestra 4 cortes distintos en Y para ver si la respuesta
% cambia a lo largo de la dimensión Y
figure('Name','Cortes en Y - Target vs Magnetización');
iy_cortes = [round(Ny*1/5), round(Ny*2/5), round(Ny*3/5), round(Ny*4/5)];

for k = 1:4
    iy = iy_cortes(k);

    subplot(4, 2, (k-1)*2 + 1);
    imagesc(x, frange, dtarget(:,:,iy), [0 flip*pi/180]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]');
    colorbar;
    title(sprintf('Target | Y = %.1f cm', y(iy)));

    subplot(4, 2, (k-1)*2 + 2);
    imagesc(x, frange, abs(mopt(:,:,iy)), [0 flip*pi/180]);
    xlabel('X [cm]'); ylabel('Frec. [Hz]');
    colorbar;
    title(sprintf('Magnetización | Y = %.1f cm', y(iy)));
end
sgtitle('Cortes en distintas posiciones Y');
colormap('parula');

%Coge 4 posiciones Y distintas (al 20%, 40%, 60% y 80% del FOV) y para cada una
% pinta el target y la magnetización en el plano frecuencia × X.
% Es como cortar el volumen 3D en rodajas horizontales. Como el target es uniforme
% en Y debería verse igual en los 4 cortes. Si la magnetización cambia entre cortes
% significa que la trayectoria espiral no es perfectamente simétrica en Y.

%% VISUALIZACIÓN 3D - Vista 2: Cortes en frecuencia (plano X-Y)
% Fijamos una frecuencia y vemos el plano espacial X-Y completo.
% Esta es la figura más importante: muestra la selectividad 2D espacial.
figure('Name','Cortes en frecuencia - plano X-Y');

freqs_interes = [fatfreq, fatfreq/2, 0, 50];
labels_freqs  = {'Banda grasa (~-220Hz)', 'Intermedia', 'Banda agua (0Hz)', 'Fuera banda'};

for k = 1:4
    [~, if_k] = min(abs(frange - freqs_interes(k)));

    subplot(4, 2, (k-1)*2 + 1);
    imagesc(x, y, squeeze(dtarget(if_k,:,:)).');
    colormap(gca,'parula'); axis image; colorbar;
    xlabel('X [cm]'); ylabel('Y [cm]');
    title(sprintf('Target | %s | f=%.0fHz', labels_freqs{k}, frange(if_k)));

    subplot(4, 2, (k-1)*2 + 2);
    imagesc(x, y, squeeze(abs(mopt(if_k,:,:))).');
    colormap(gca,'parula'); axis image; colorbar;
    xlabel('X [cm]'); ylabel('Y [cm]');
    title(sprintf('Magnetización | %s', labels_freqs{k}));
end
sgtitle('Plano X-Y en distintas frecuencias');

%Ahora al revés: fijas una frecuencia y ves el plano espacial X × Y.
% Se Elige 4 frecuencias interesantes: la banda de grasa (~-220 Hz), una frecuencia
% intermedia, la banda de agua (0 Hz) y una fuera de banda (50 Hz).
%squeeze elimina la dimensión de tamaño 1 que queda al coger una sola frecuencia,
% para que quede un array 2D de (Nx × Ny). La trasposición .' es para que X quede
% en el eje horizontal e Y en el vertical.
%En la banda de grasa se debería ver franjas verticales en el target (las columnas dfat)
% y la magnetización intentando reproducirlas. En la banda de agua todo debería
% ser cero porque no se quiere excitar el agua.

%% VISUALIZACIÓN 3D - Vista 3: Perfil espectral en puntos X-Y
% Fijamos un punto (X,Y) y vemos su respuesta a cada frecuencia.
% Comparamos un punto de grasa con uno de agua.
figure('Name','Perfil espectral en puntos X-Y');

% Punto en zona de grasa (columna 15 está en dfat=[9:19])
ix_fat = 15;  iy_fat = round(Ny/2);

% Punto en zona de agua (columna 6 está en dwater=[4:8])
ix_wat = 6;   iy_wat = round(Ny/2);

plot(frange, abs(mopt(:, ix_fat, iy_fat)), 'r-',  'LineWidth', 2, ...
    'DisplayName', sprintf('Magnetización grasa  X=%.1f Y=%.1f cm', x(ix_fat), y(iy_fat)));
hold on;
plot(frange, abs(mopt(:, ix_wat, iy_wat)), 'b-',  'LineWidth', 2, ...
    'DisplayName', sprintf('Magnetización agua   X=%.1f Y=%.1f cm', x(ix_wat), y(iy_wat)));
plot(frange, dtarget(:, ix_fat, iy_fat),   'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'Target grasa');
plot(frange, dtarget(:, ix_wat, iy_wat),   'b--', 'LineWidth', 1.5, ...
    'DisplayName', 'Target agua');
xline(fatfreq, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Frecuencia grasa');
xline(0,       'g:', 'LineWidth', 1.5, 'DisplayName', 'Frecuencia agua');
hold off;
xlabel('Frecuencia [Hz]'); ylabel('|M_{xy}|');
title('Perfil espectral: punto de grasa vs punto de agua');
legend('Location','northwest'); grid on;
ylim([0 flip*pi/180 * 1.5]);

%Se Fija un punto concreto del espacio (una columna X y una fila Y) y se ve cómo
% responde a cada frecuencia del rango. Es como preguntar "en este punto del espacio,
% ¿qué frecuencias excita el pulso?".
%Se comparan dos puntos: uno en zona de grasa (columna 15, dentro de dfat=[9:19])
% y otro en zona de agua (columna 6, dentro de dwater=[4:8]).
% En el punto de grasa se debería ver un pico alto en -220 Hz y bajo en 0 Hz.
% En el punto de aguao  debería ser bajo en todas las frecuencias.
% Si los dos puntos se comportan de forma muy parecida, significa que SPSP
% no está consiguiendo distinguir bien la posición espacial.

%Resultado: Que la gráfica muestra que con los parámetros actuales (Nx=Ny=32, flip=30°)
% el pulso SPSP no consigue distinguir bien entre la posición de grasa y la de agua.
% El target sí está bien definido (pico rojo discontinuo), pero la magnetización
% conseguida es demasiado baja y uniforme. Esto es una limitación del diseño del pulso
% con resolución espacial reducida, no un error del código.