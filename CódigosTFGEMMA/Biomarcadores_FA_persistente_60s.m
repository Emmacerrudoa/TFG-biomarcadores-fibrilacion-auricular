clear
clc
close all

%% ============================================================
% EXTRACCION DE BIOMARCADORES: FA PERSISTENTE
%
% Se parte de las ventanas .mat de 120 s ya construidas durante el
% preprocesamiento y libres de artefactos.
%
% Para cada paciente se selecciona aleatoriamente una pareja valida
% formada por dos ventanas consecutivas de 120 s.
%
% Con ellas se forma un tramo continuo de 240 s y se analizan:
%
%   MOMENTO 1:   0 a  60 s
%   DESCANSO :  60 a 120 s
%   MOMENTO 2: 120 a 180 s
%
% Cada paciente aporta dos filas al Excel.
%
% Biomarcadores por ventana:
%   RR_mean, SDNN, RMSSD, SDSD, pNN20, pNN50, CV_RR,
%   SD1, SD2, SD1_SD2,
%   DF_completo_Hz y DF_residual_Hz.
%
% No se calculan LF, HF, LF/HF, LFnu, HFnu, SampEn,
% ni biomarcadores morfologicos de ondas P y T.
%% ============================================================

%% ============================
% CONFIGURACION
%% ============================

% Semilla fija para que la seleccion aleatoria sea reproducible.
rng(42)

% Modificar estas rutas segun la ubicacion local de las carpetas.
carpeta_entrada = ...
    'C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\FA_PERSISTENTE';

carpeta_out = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_60s_FA_PERSISTENTE';

if ~exist(carpeta_out,'dir')
    mkdir(carpeta_out)
end

archivo_salida = fullfile( ...
    carpeta_out, ...
    'biomarcadores_FA_persistente_60s.xlsx');

Fs_objetivo = 500;

% Ventanas dentro del tramo continuo formado por dos ventanas de 120 s.
ventanas_rel = [ ...
      0   60; ...
    120  180];

etiquetas_momento = [ ...
    "FA_persistente_minuto_1"; ...
    "FA_persistente_minuto_3"];

orden_momento = [1; 2];

% Controles minimos de calidad.
min_R_60s = 20;
min_RR_validos = 20;

%% ============================
% LOCALIZAR ARCHIVOS .MAT
%% ============================

if ~exist(carpeta_entrada,'dir')
    error('No existe la carpeta de entrada:\n%s',carpeta_entrada);
end

archivos = dir(fullfile(carpeta_entrada,'**','*.mat'));

if isempty(archivos)
    error('No se encontraron archivos .mat en:\n%s',carpeta_entrada);
end

fprintf('\nArchivos .mat encontrados: %d\n',numel(archivos));

%% ============================
% LEER METADATOS DE TODAS LAS VENTANAS
%% ============================

Registros = table;

for i = 1:numel(archivos)

    ruta = fullfile(archivos(i).folder,archivos(i).name);

    try
        S = load(ruta);

        [ecg,Fs,locs_R,id_paciente,t_ini,t_fin] = ...
            leer_ventana_guardada(S);

        if isempty(ecg) || isempty(Fs) || isempty(locs_R)
            continue
        end

        if abs(Fs-Fs_objetivo) > 1e-6
            continue
        end

        nueva = table;
        nueva.Ruta = string(ruta);
        nueva.Paciente = string(id_paciente);
        nueva.t_ini = t_ini;
        nueva.t_fin = t_fin;
        nueva.Duracion = numel(ecg)/Fs;

        Registros = [Registros; nueva]; %#ok<AGROW>

    catch ME
        fprintf('No se pudo leer %s: %s\n',ruta,ME.message);
    end
end

if isempty(Registros)
    error('No se pudo leer ninguna ventana valida.');
end

Registros = sortrows(Registros,{'Paciente','t_ini'});

fprintf('Ventanas validas leidas: %d\n',height(Registros));
fprintf('Pacientes distintos: %d\n',numel(unique(Registros.Paciente)));

%% ============================
% SELECCIONAR ALEATORIAMENTE DOS
% VENTANAS CONSECUTIVAS POR PACIENTE
%% ============================

pacientes = unique(Registros.Paciente,'stable');

filas_resultados = {};
pacientes_conservados = strings(0,1);

for p = 1:numel(pacientes)

    paciente = pacientes(p);

    idx_p = Registros.Paciente == paciente;
    Rp = Registros(idx_p,:);
    Rp = sortrows(Rp,'t_ini');

    pareja_encontrada = false;

    %% BUSCAR TODAS LAS PAREJAS TEMPORALMENTE CONSECUTIVAS

    tolerancia = 1/Fs_objetivo;
    indices_parejas = [];

    for k = 1:height(Rp)-1

        if abs(Rp.t_ini(k+1)-Rp.t_fin(k)) <= tolerancia
            indices_parejas(end+1,1) = k; %#ok<AGROW>
        end
    end

    if isempty(indices_parejas)

        fprintf( ...
            'Paciente %s sin parejas de ventanas consecutivas.\n', ...
            paciente);

        continue
    end

    %% CAMBIAR ALEATORIAMENTE EL ORDEN DE LAS PAREJAS

    orden_aleatorio = randperm(numel(indices_parejas));
    indices_parejas = indices_parejas(orden_aleatorio);

    %% PROBAR LAS PAREJAS EN ORDEN ALEATORIO
    % Se conserva la primera que permite calcular los dos momentos.

    for pp = 1:numel(indices_parejas)

        k = indices_parejas(pp);

        S1 = load(Rp.Ruta(k));
        S2 = load(Rp.Ruta(k+1));

        [ecg1,Fs1,R1,id1,t_ini1,t_fin1] = ...
            leer_ventana_guardada(S1);

        [ecg2,Fs2,R2,id2,t_ini2,t_fin2] = ...
            leer_ventana_guardada(S2);

        if abs(Fs1-Fs_objetivo) > 1e-6 || ...
                abs(Fs2-Fs_objetivo) > 1e-6
            continue
        end

        if string(id1) ~= string(id2)
            continue
        end

        % Unir ECG.
        ecg240 = [ecg1(:); ecg2(:)];

        % Reubicar los R de la segunda ventana.
        R1 = round(R1(:));
        R2 = round(R2(:)) + numel(ecg1);

        R240 = unique([R1;R2]);

        R240 = R240( ...
            R240 >= 1 & ...
            R240 <= numel(ecg240));

        % Comprobar duracion suficiente.
        if numel(ecg240) < round(180*Fs_objetivo)
            continue
        end

        resultados_paciente = cell(2,1);
        paciente_completo = true;

        %% ANALIZAR LOS DOS MOMENTOS

        for mm = 1:2

            ini60 = round(ventanas_rel(mm,1)*Fs_objetivo)+1;
            fin60 = round(ventanas_rel(mm,2)*Fs_objetivo);

            ecg60 = ecg240(ini60:fin60);
            ecg60 = ecg60(:);

            R60 = R240( ...
                R240 >= ini60 & ...
                R240 <= fin60)-ini60+1;

            R60 = unique(round(R60(:)));

            if numel(R60) < min_R_60s
                paciente_completo = false;
                break
            end

            B = calcular_biomarcadores_FA_60s( ...
                ecg60, ...
                R60, ...
                Fs_objetivo, ...
                min_RR_validos);

            if ~B.Valida
                paciente_completo = false;
                break
            end

            resultados_paciente{mm} = B;
        end

        if ~paciente_completo
            continue
        end

        pacientes_conservados(end+1,1) = paciente; %#ok<SAGROW>

        %% GUARDAR DOS FILAS POR PACIENTE

        for mm = 1:2

            B = resultados_paciente{mm};

            filas_resultados(end+1,:) = { ... %#ok<SAGROW>
                paciente, ...
                t_ini1, ...
                t_fin2, ...
                etiquetas_momento(mm), ...
                orden_momento(mm), ...
                ventanas_rel(mm,1), ...
                ventanas_rel(mm,2), ...
                string(Rp.Ruta(k)), ...
                string(Rp.Ruta(k+1)), ...
                B.N_R, ...
                B.N_RR, ...
                B.RR_mean, ...
                B.SDNN, ...
                B.RMSSD, ...
                B.SDSD, ...
                B.pNN20, ...
                B.pNN50, ...
                B.CV_RR, ...
                B.SD1, ...
                B.SD2, ...
                B.SD1_SD2, ...
                B.DF_completo_Hz, ...
                B.DF_residual_Hz};
        end

        fprintf( ...
            ['Paciente %s | pareja seleccionada aleatoriamente: ' ...
             '%.1f a %.1f s\n'], ...
            paciente,t_ini1,t_fin2);

        pareja_encontrada = true;
        break
    end

    if ~pareja_encontrada

        fprintf( ...
            ['Paciente %s sin ninguna pareja aleatoria que permita ' ...
             'calcular los dos momentos.\n'], ...
            paciente);
    end
end

%% ============================
% CREAR Y GUARDAR TABLA
%% ============================

if isempty(filas_resultados)
    error(['No se obtuvo ningun paciente con FA persistente ' ...
        'y dos ventanas validas.']);
end

nombres_columnas = { ...
    'Paciente', ...
    'Tiempo_inicio_bloque_s', ...
    'Tiempo_fin_bloque_s', ...
    'Momento', ...
    'Orden_momento', ...
    'Ventana_ini_rel_s', ...
    'Ventana_fin_rel_s', ...
    'Archivo_ventana_1', ...
    'Archivo_ventana_2', ...
    'N_R', ...
    'N_RR', ...
    'RR_mean', ...
    'SDNN', ...
    'RMSSD', ...
    'SDSD', ...
    'pNN20', ...
    'pNN50', ...
    'CV_RR', ...
    'SD1', ...
    'SD2', ...
    'SD1_SD2', ...
    'DF_completo_Hz', ...
    'DF_residual_Hz'};

T = cell2table( ...
    filas_resultados, ...
    'VariableNames',nombres_columnas);

T = sortrows( ...
    T, ...
    {'Paciente','Orden_momento'});

writetable(T,archivo_salida);

fprintf('\nTabla guardada en:\n%s\n',archivo_salida);
fprintf('Pacientes conservados: %d\n',numel(unique(T.Paciente)));
fprintf('Filas totales: %d\n',height(T));
fprintf('Cada paciente aporta dos filas: minuto 1 y minuto 3.\n');
fprintf('La pareja se ha seleccionado aleatoriamente con rng(42).\n');

%% ============================================================
% FUNCIONES LOCALES
%% ============================================================

function [ecg,Fs,locs_R,id_paciente,t_ini,t_fin] = ...
    leer_ventana_guardada(S)

ecg = double(S.ventana);
ecg = ecg(:);

Fs = double(S.Fs);

locs_R = double(S.locs_R);
locs_R = locs_R(:);

id_paciente = string(S.nombre_registro);

t_ini = double(S.t_ini);
t_fin = double(S.t_fin);

end

function B = calcular_biomarcadores_FA_60s( ...
    ecg,R,Fs,min_RR)

B = struct( ...
    'N_R',numel(R), ...
    'N_RR',0, ...
    'RR_mean',NaN, ...
    'SDNN',NaN, ...
    'RMSSD',NaN, ...
    'SDSD',NaN, ...
    'pNN20',NaN, ...
    'pNN50',NaN, ...
    'CV_RR',NaN, ...
    'SD1',NaN, ...
    'SD2',NaN, ...
    'SD1_SD2',NaN, ...
    'DF_completo_Hz',NaN, ...
    'DF_residual_Hz',NaN, ...
    'Valida',false);

%% BIOMARCADORES RR

RR = diff(R)/Fs;

RR = RR(isfinite(RR));
RR = RR(RR >= 0.30 & RR <= 1.50);

B.N_RR = numel(RR);

if numel(RR) < min_RR
    return
end

B.RR_mean = mean(RR,'omitnan');
B.SDNN = std(RR,0,'omitnan');

if B.RR_mean > 0
    B.CV_RR = B.SDNN/B.RR_mean;
end

dRR = diff(RR);

if isempty(dRR)
    return
end

B.RMSSD = sqrt(mean(dRR.^2,'omitnan'));
B.SDSD = std(dRR,0,'omitnan');

B.pNN20 = 100*mean(abs(dRR) > 0.020);
B.pNN50 = 100*mean(abs(dRR) > 0.050);

%% DIAGRAMA DE POINCARE

RRn = RR(1:end-1);
RRn1 = RR(2:end);

xp = (RRn + RRn1)/sqrt(2);
yp = (RRn1 - RRn)/sqrt(2);

B.SD1 = std(yp,0,'omitnan');
B.SD2 = std(xp,0,'omitnan');

if isfinite(B.SD2) && B.SD2 > 0
    B.SD1_SD2 = B.SD1/B.SD2;
end

%% FRECUENCIA DOMINANTE EN FA
% La senal almacenada ya procede del preprocesamiento general 0.5-40 Hz.

[B.DF_completo_Hz,~,~,~] = ...
    frecuencia_dominante2_FA120s(ecg,Fs);

residual = cancelar_QRST_plantilla_medianaFA( ...
    ecg,R,Fs);

if ~isempty(residual) && all(isfinite(residual))

    [B.DF_residual_Hz,~,~,~] = ...
        frecuencia_dominante2_FA120s(residual,Fs);
end

%% CONTROL FINAL

valores = [ ...
    B.RR_mean, ...
    B.SDNN, ...
    B.RMSSD, ...
    B.SDSD, ...
    B.pNN20, ...
    B.pNN50, ...
    B.CV_RR, ...
    B.SD1, ...
    B.SD2, ...
    B.SD1_SD2, ...
    B.DF_completo_Hz, ...
    B.DF_residual_Hz];

B.Valida = all(isfinite(valores));
end

function [DF, f_axis, Pxx, peak_power] = ...
    frecuencia_dominante2_FA120s(x, Fs)

DF = NaN;
peak_power = NaN;
f_axis = [];
Pxx = [];

x = x(:);

if isempty(Fs) || ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
    return
end

win = round(30 * Fs);

if isempty(x) || numel(x) < win
    return
end

if any(~isfinite(x))
    return
end

%% Quitar media y tendencia

x = x - mean(x, 'omitnan');
x = detrend(x);

%% Espectro de potencia con Welch

% Segmentos de 30 s
% Solapamiento del 50 %

noverlap = round(0.5 * win);
Nfft = 2^nextpow2(win);

[Pxx_raw, f_axis] = pwelch( ...
    x, ...
    hamming(win), ...
    noverlap, ...
    Nfft, ...
    Fs);

%% Suavizar el espectro completo

Pxx = movmean(Pxx_raw, 5);

%% Buscar pico dominante entre 3 y 9 Hz

idx = f_axis >= 3 & f_axis <= 9;

if ~any(idx)
    return
end

P_banda = Pxx(idx);
f_valid = f_axis(idx);

[Pmax, im] = max(P_banda);

DF = f_valid(im);
peak_power = P_banda(im);

%% Control de calidad

media_banda = mean(P_banda, 'omitnan');

if ~isfinite(Pmax) || ~isfinite(media_banda) || media_banda <= 0
    DF = NaN;
    peak_power = NaN;
    return
end

if Pmax < 1.5 * media_banda
    DF = NaN;
    peak_power = NaN;
    return
end

%% Rechazar picos proximos a los bordes

if DF <= 3.1 || DF >= 8.9
    DF = NaN;
    peak_power = NaN;
    return
end

end

function [x_residual, plantilla, t_plantilla, latidos_validos] = cancelar_QRST_plantilla_medianaFA(x, locs_R, Fs)

% CANCELAR_QRST_PLANTILLA_MEDIANAFA
% Cancela los complejos QRS-T mediante sustraccion de una plantilla mediana.
%
% Version para FA:
%   - 200 ms antes del R
%   - 450 ms despues del R
%
% Entrada:
%   x      -> senal ECG filtrada
%   locs_R -> posiciones de los picos R en muestras
%   Fs     -> frecuencia de muestreo
%
% Salida:
%   x_residual      -> senal residual con menor contribucion QRS-T
%   plantilla       -> plantilla mediana QRST base
%   t_plantilla     -> eje temporal de la plantilla respecto al pico R, en segundos
%   latidos_validos -> latidos empleados para construir la plantilla

x = x(:);
x_residual = x;

plantilla = [];
t_plantilla = [];
latidos_validos = [];

if isempty(x) || isempty(locs_R)
    x_residual = [];
    return
end

locs_R = limpiar_locs_local(locs_R, length(x));

if numel(locs_R) < 5
    x_residual = [];
    return
end

%% 1) DEFINIR VENTANA QRS-T ALREDEDOR DE CADA R

pre_R  = round(0.20 * Fs);   % 200 ms antes del R
post_R = round(0.45 * Fs);   % 450 ms despues del R

L = pre_R + post_R + 1;
t_plantilla = (-pre_R:post_R) / Fs;

latidos = nan(numel(locs_R), L);
validos = false(numel(locs_R), 1);

%% 2) EXTRAER LATIDOS ALINEADOS POR R

for i = 1:numel(locs_R)

    R = locs_R(i);

    ini = R - pre_R;
    fin = R + post_R;

    if ini < 1 || fin > length(x)
        continue
    end

    latido = x(ini:fin);

    if any(~isfinite(latido))
        continue
    end

    latidos(i,:) = latido(:)';
    validos(i) = true;

end

latidos_validos = latidos(validos,:);

if size(latidos_validos,1) < 5
    x_residual = [];
    plantilla = [];
    t_plantilla = [];
    latidos_validos = [];
    return
end

%% 3) CREAR PLANTILLA MEDIANA QRS-T

plantilla = median(latidos_validos, 1, 'omitnan');

%% 4) RESTAR LA PLANTILLA EN CADA LATIDO

for i = 1:numel(locs_R)

    R = locs_R(i);

    ini = R - pre_R;
    fin = R + post_R;

    if ini < 1 || fin > length(x_residual)
        continue
    end

    segmento = x_residual(ini:fin);

    if any(~isfinite(segmento))
        continue
    end

    % Ajuste de amplitud para adaptar la plantilla a cada latido.
    num = segmento(:)' * plantilla(:);
    den = plantilla(:)' * plantilla(:);

    if den > 0
        escala = num / den;
    else
        escala = 1;
    end

    plantilla_ajustada = escala * plantilla(:);

    x_residual(ini:fin) = segmento(:) - plantilla_ajustada;

end

%% 5) CENTRAR SENAL RESIDUAL

x_residual = x_residual - mean(x_residual, 'omitnan');

end

function locs = limpiar_locs_local(locs, N)

if isempty(locs)
    locs = [];
    return
end

locs = round(locs(:));
locs = locs(isfinite(locs));
locs = unique(locs);
locs = locs(locs >= 1 & locs <= N);

end