% stir_phantom.m
% Ajustado al phantom con grasa en capa 3 (anillo exterior)

%% BLOQUE 0 - Cargar y visualizar
load('/Users/sonita/Desktop/4º/TFG/Datos TFG/phantomdigital.mat');  % carga phantom (128x128) y tejidos (128x128x11)

figure('Name', 'Phantom original');
imagesc(phantom); colormap(gca,'parula'); axis image; colorbar;
title('Phantom digital');
xlabel('X [píxeles]'); ylabel('Y [píxeles]');

%% BLOQUE 1 - Parámetros del sistema y T1/T2
Sys = sys_sola();
T1  = t1_values_1_5;
T2  = t2_values_1_5;

fprintf('\n=== T1 de tejidos a 1.5T ===\n');
fprintf('Grasa:       %6.1f ms  <- T1 CORTO, STIR lo suprime\n', T1.fat);
fprintf('Músculo:     %6.1f ms\n', T1.muscle);
fprintf('Agua/CSF:    %6.1f ms  <- T1 LARGO\n', T1.csf);
fprintf('Hígado:      %6.1f ms\n', T1.liver);
fprintf('Mat. blanca: %6.1f ms\n', T1.wm);
fprintf('Mat. gris:   %6.1f ms\n', T1.gm);

% Asignación de T1 a cada capa según lo que se ve en la figura de tejidos:
% Tejido 1  = fondo/aire         -> T1 = 0 (no hay señal)
% Tejido 2  = interior del cuerpo -> músculo/tejido blando
% Tejido 3  = anillo exterior    -> GRASA
% Tejido 4  = óvalo derecho      -> músculo
% Tejido 5  = óvalo izquierdo    -> músculo
% Tejido 6  = círculo superior   -> CSF/líquido
% Tejidos 7-11 = puntos pequeños -> sangre

T1_por_capa = [
    0,           ...  % capa  1: fondo/aire (sin señal)
    T1.muscle,   ...  % capa  2: interior general del cuerpo
    T1.fat,      ...  % capa  3: GRASA (anillo exterior) <- clave STIR
    T1.muscle,   ...  % capa  4: óvalo derecho
    T1.muscle,   ...  % capa  5: óvalo izquierdo
    T1.csf,      ...  % capa  6: círculo superior (líquido)
    T1.muscle,   ...  % capa  7: punto pequeño 
    T1.muscle,   ...  % capa  8: punto pequeño
    T1.muscle,   ...  % capa  9: punto pequeño
    T1.muscle,   ...  % capa 10: punto pequeño
    T1.muscle    ...  % capa 11: punto pequeño
];  % [ms]

% Construir mapa T1 píxel a píxel
T1map = zeros(128,128);
for k = 1:11
    T1map = T1map + tejidos(:,:,k) * T1_por_capa(k);
end

%% Mapa T1 con colores más visibles
figure('Name','Mapa T1');

% Para que el fondo no domine la escala de color
T1map_plot = T1map;
T1map_plot(T1map_plot == 0) = NaN;

hImg = imagesc(T1map_plot);
set(hImg, 'AlphaData', ~isnan(T1map_plot));

axis image;
set(gca,'Color',[0.85 0.85 0.85]);   % fondo gris claro

% Escala de color solo con tejidos reales, no con el fondo
valsT1 = T1map(T1map > 0);
caxis([min(valsT1) max(valsT1)]);

colormap(gca,'parula');
h = colorbar;
ylabel(h,'T_1 [ms]');

hold on;
contour(double(tejidos(:,:,3)), [0.5 0.5], 'w-', 'LineWidth', 1.5);
hold off;

title('Mapa T_1 [ms] del phantom');
xlabel('X [px]');
ylabel('Y [px]');

%% BLOQUE 2 - Máscara de grasa (capa 3)
capa_grasa    = 3;
mascara_grasa = logical(tejidos(:,:, capa_grasa));

% Verificar visualmente que la máscara coincide con el anillo amarillo
figure('Name', 'Verificación máscara de grasa', ...
       'Position', [100 100 1000 420]);

subplot(1,2,1);
imagesc(phantom); 
colormap(gca,'parula'); 
axis image; 
colorbar;
title('Phantom original');

% Guardar posición real de la primera imagen después del colorbar
ax1 = gca;
pos1 = ax1.Position;

subplot(1,2,2);
phantom_norm = phantom / max(phantom(:));
img_rgb      = repmat(phantom_norm, [1 1 3]);

img_rgb(:,:,1) = min(img_rgb(:,:,1) + 0.7 * double(mascara_grasa), 1);
img_rgb(:,:,2) = img_rgb(:,:,2) .* (1 - double(mascara_grasa));
img_rgb(:,:,3) = img_rgb(:,:,3) .* (1 - double(mascara_grasa));

imshow(img_rgb); 
axis image;
title({'Grasa marcada en rojo', 'capa 3 / anillo exterior'});

% Hacer la imagen de la derecha un poco más pequeña
ax2 = gca;
pos2 = ax2.Position;

escala = 0.90; 

nuevo_ancho = pos1(3) * escala;
nuevo_alto  = pos1(4) * escala;

ax2.Position = [pos2(1), ...
                pos1(2) + (pos1(4)-nuevo_alto)/2, ...
                nuevo_ancho, ...
                nuevo_alto];

fprintf('\nPíxeles de grasa: %d (%.1f%% del total)\n', ...
    sum(mascara_grasa(:)), 100*mean(mascara_grasa(:)));

%% BLOQUE 3 - Curvas de recuperación T1 y TI óptimo
T1_fat_s = T1.fat    / 1000;   % [s]
T1_mus_s = T1.muscle / 1000;
T1_csf_s = T1.csf    / 1000;

% TI óptimo: momento en que Mz de grasa = 0
% De:  0 = 1 - 2*exp(-TI/T1)  =>  TI = T1 * ln(2)
TI_fat    = T1_fat_s * log(2);
TI_muscle = T1_mus_s * log(2);
TI_water  = T1_csf_s * log(2);

fprintf('\n=== Tiempos de inversión TI = T1 * ln(2) ===\n');
fprintf('TI grasa:   %6.1f ms  <- usaremos ESTE para STIR\n', TI_fat*1000);
fprintf('TI músculo: %6.1f ms\n', TI_muscle*1000);
fprintf('TI agua:    %6.1f ms\n', TI_water*1000);

% Curvas de recuperación para ver dónde cada tejido cruza por cero
time_rec      = linspace(0, 3, 3000);   % [s]
Mz_fat_rec    = 1 - 2*exp(-time_rec / T1_fat_s);
Mz_muscle_rec = 1 - 2*exp(-time_rec / T1_mus_s);
Mz_water_rec  = 1 - 2*exp(-time_rec / T1_csf_s);

figure('Name', 'Recuperación T1 - TI óptimo para grasa');

plot(time_rec*1000, Mz_fat_rec, 'r-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Grasa (T_1 = %.0f ms)', T1.fat));
hold on;

plot(time_rec*1000, Mz_muscle_rec, 'g-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Músculo (T_1 = %.0f ms)', T1.muscle));

plot(time_rec*1000, Mz_water_rec, 'b-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Agua/CSF (T_1 = %.0f ms)', T1.csf));

% Línea horizontal donde Mz = 0
yline(0, 'k--', 'LineWidth', 1.2, ...
    'DisplayName', 'M_z = 0');

% Línea vertical del TI de la grasa
xline(TI_fat*1000, 'r--', ...
    sprintf('TI_{grasa} = %.0f ms', TI_fat*1000), ...
    'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', ...
    'LabelHorizontalAlignment', 'left', ...
    'DisplayName', sprintf('TI grasa = %.0f ms', TI_fat*1000));

hold off;

xlabel('Tiempo [ms]');
ylabel('M_z');
title('Recuperación T_1 tras pulso de inversión 180°');
legend('Location','southeast');
grid on;
ylim([-1.1 1.1]);
% En TI_fat la curva roja cruza cero -> grasa no da señal en ese instante

%% BLOQUE 4 - Aplicar STIR al phantom
% Señal STIR de cada píxel:
%   S(x,y) = |Mz(TI)| * phantom(x,y)
%   Mz(TI) = 1 - 2*exp(-TI / T1(x,y))
% Con TI = TI_fat:
%   grasa   -> Mz ~ 0  -> señal suprimida
%   músculo -> Mz > 0  -> señal conservada 
%   agua    -> Mz < 0  -> señal conservada 

T1map_s   = T1map / 1000;           % convertir a segundos
Mz_en_TI  = zeros(128,128);
mask_valido = T1map_s > 0;          % excluir fondo donde T1 = 0

Mz_en_TI(mask_valido) = 1 - 2*exp(-TI_fat ./ T1map_s(mask_valido));

imagen_STIR   = abs(Mz_en_TI) .* phantom;
imagen_normal = phantom;

% Comprobar que la grasa queda en cero
Mz_grasa_en_TI = 1 - 2*exp(-TI_fat / T1_fat_s);
fprintf('\nMz de grasa en TI = %.0f ms: %.6f (ideal = 0)\n', ...
    TI_fat*1000, Mz_grasa_en_TI);

%% BLOQUE 5 - Figura principal de resultados
figure('Name','Resultado STIR', 'Position', [50 50 1400 420]);

subplot(1,4,1);
imagesc(imagen_normal); colormap(gca,'gray'); axis image; colorbar;
title('Sin STIR (imagen normal)');
xlabel('X [px]'); ylabel('Y [px]');

subplot(1,4,2);
imagesc(imagen_STIR); colormap(gca,'gray'); axis image; colorbar;
title(sprintf('Con STIR  (TI = %.0f ms)', TI_fat*1000));
xlabel('X [px]');
% El anillo exterior debe aparecer negro (suprimido)

subplot(1,4,3);
imagesc(imagen_normal - imagen_STIR); colormap(gca,'hot'); axis image; colorbar;
title('Señal suprimida por STIR');
xlabel('X [px]');
% Debe verse brillante SOLO en el anillo de grasa

subplot(1,4,4);
imagesc(imagen_STIR); colormap(gca,'gray'); axis image; colorbar;
hold on;
% Dibujar contorno de la máscara de grasa en rojo para verificar
contour(double(mascara_grasa), 1, 'r-', 'LineWidth', 1.5);
hold off;
title('STIR + contorno grasa (rojo)');
xlabel('X [px]');

%% BLOQUE 6 - Cuantificación numérica
senal_grasa_normal = mean(imagen_normal(mascara_grasa));
senal_grasa_STIR   = mean(imagen_STIR(mascara_grasa));
supresion_pct      = 100 * (1 - senal_grasa_STIR / senal_grasa_normal);

% Señal en músculo (capa 2, interior general)
mascara_musculo      = logical(tejidos(:,:,2));
senal_musculo_normal = mean(imagen_normal(mascara_musculo));
senal_musculo_STIR   = mean(imagen_STIR(mascara_musculo));

fprintf('\n=== Resultado cuantitativo ===\n');
fprintf('%-25s %-12s %-12s %-15s\n', 'Tejido', 'Sin STIR', 'Con STIR', '% supresión');
fprintf('%-25s %-12.4f %-12.4f %-14.1f%%\n', ...
    'Grasa (capa 3)', senal_grasa_normal, senal_grasa_STIR, supresion_pct);
fprintf('%-25s %-12.4f %-12.4f %-14.1f%%\n', ...
    'Músculo (capa 2)', senal_musculo_normal, senal_musculo_STIR, ...
    100*(1 - senal_musculo_STIR/senal_musculo_normal));
fprintf('\nMz grasa en TI_fat: %.6f (ideal = 0 -> 100%% supresión)\n', ...
    Mz_grasa_en_TI);