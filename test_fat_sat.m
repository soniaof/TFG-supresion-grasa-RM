% test_fat_sat.m
% añadir directorios
codedir=cd;                 % directorio donde está el código
datdir='C:\Users\sydney.williams\OneDrive - Universidad Rey Juan Carlos\Documentos\Students\TFG\2026\SoniaODLF\Datos\';
dirbloch='C:\Users\sydney.williams\OneDrive - Universidad Rey Juan Carlos\Documentos\Students\TFG\2026\SoniaODLF\Codigo\BlochSim\'; % directorio donde está el simulador de Bloch
addpath(genpath(dirbloch));
dirirt='C:\Users\sydney.williams\OneDrive - Universidad Rey Juan Carlos\Documentos\ResearchCode\fessler\irt\';                      % directorio donde está el toolbox IRT
cd(dirirt); setup; cd(codedir);
%% Crear parámetros para simulacion 

% variables que pertenecen at 1.5T
Sys=sys_sola;                               % carga parámetros del sistema, Siemens Sola Fit 1.5T
T1=t1_values_1_5;                           % carga valores del T1 a 1.5T
T2=t2_values_1_5;                           % carga valores del T2 a 1.5T
% pulso RF
Trf=5e-3;                                   % tiempo del pulso RF [sec]
time=[0:Sys.dt:Trf];                        % puntos del tiempo de pulse RF [sec]
Npts=length(time);                          % longitud de pulso
flip=90;                                    % angulo de giro [deg]
b1amp=(flip*pi/180)/(2*pi*Sys.gamma*Trf);   % amplitud de pulso de RF [***calcular unidades***]
b1rf=[0; b1amp*ones((Npts-2),1); 0];        % crear pulso RF "hard" (rectangular)
gx=zeros(Npts,1); gy=gx; gz=gx;               % sin gradientes para pulso "hard"

%% Representar pulsos y gradientes
figure; 
subplot(2,1,1);
yyaxis left
plot(time,abs(b1rf)); axis([time(1) time(end) 0 round(max(abs(b1rf))*1e7)/1e7]);
ylabel('|B1| (**unidades**)')
yyaxis right
plot(time,angle(b1rf)); axis([time(1) time(end) -pi pi]);
xlabel('Time (sec)'); ylabel('<B1 (rad)');
subplot(2,1,2);
plot(time,gx); hold on; plot(time,gy); plot(time,gz);
hold off; xlabel('Time (sec)'); ylabel('G_{x,y,z} (mT/m)');

%% Simulador de Bloch
FOVx=25;                    % [cm]
Nx=128;
x=[-Nx/2:Nx/2-1]*FOVx/Nx; x=x(:);
freq=zeros(size(x));        % [Hz]
sens=ones([1 size(x)]);     % senisbilidad de RF, no lo utilizamos
% ***ojo: para correr el simulador de Bloch, hay que convertir las unidades 1T = 1e4 G ***
% b1rfsim = ?
% gxsim = ? gysim = ? gzsim = ?
% T1sim = ? T2sim = ? (elegir un tejido, por ejemplo grasa)
[Mx,My,Mz] = blochCim(b1rfsim,[gxsim, gysim, gzsim],Sys.dt,T1sim,T2sim,freq,x,0,sens);
Mxy=Mx+1i*My;

%% Representar simulación de Bloch "on-resonance" (a la frecuencia de agua)
figure;
subplot(2,1,1);
yyaxis left;
plot(x,abs(Mxy)); axis([x(1) x(end) -1 1]);
ylabel('|M_{xy}|'); title('M_{xy}');
yyaxis left;
plot(x,angle(Mxy)); axis([x(1) x(end) -pi pi]);
xlabel('Time (sec)'); ylabel('<M_{xy}');
subplot(2,1,2);
plot(x,Mz); axis([x(1) x(end) -1 1]);
xlabel('Time (sec)'); ylabel('|M_{xy}|'); title('M_{z}');

%% **repetir la simulación del pulso hard con la frecuencia de grasa (calcularla desde el ppm en sys_sola.m)**

% freqf = ?
% [Mxf,Myf,Mzf] = blochCim(...

%% **cambiar el pulso para saturar la frecuencia de grasa y simular **
% b1fat=b1rf.*exp(2*pi*1i*time(:)*freqf); 
% [Mxfsat,Myfsat,Mzfsat] = blochCim(...

%% **crear una simulación multi-epsectral para ver la saturación de frecuencias**

frange=[-300:10:100];                                   % rango de frecuencias [Hz]
Nf=length(frange);                                      % número de frecuencias
Mxspec=zeros(Nx,Nf); Myspec=Mxspec; Mzspec=Mxspec;      % Marcador de posición
Mxyspec=Mxspec;
% for ii=1:Nf
%     [Mxspec(:,ii),Myspec(:,ii),Mzspec(:,ii)] = blochCim(...
%     Mxyspec(:,ii)=Mxspec(:,ii)+1i*Myspec(:,ii);
% end

%% **crear una simulación temporal para ver la evolución de magnetización con la aplicación de RF**