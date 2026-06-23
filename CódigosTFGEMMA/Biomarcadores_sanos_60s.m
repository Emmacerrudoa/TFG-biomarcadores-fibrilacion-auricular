clear
clc
close all

%% ============================================================
% EXTRACCION DE BIOMARCADORES: SUJETOS SANOS
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
%   DF_completo_Hz y DF_residual_Hz,
%   biomarcadores morfologicos de las ondas P y T.
%
% No se calculan LF, HF, LF/HF, LFnu, HFnu ni SampEn.
%% ============================================================

%% ============================
% CONFIGURACION
%% ============================

% Semilla fija para que la seleccion aleatoria sea reproducible.
rng(42)

% Modificar estas rutas segun la ubicacion local de las carpetas.
carpeta_entrada = ...
    'C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\SANO';

carpeta_out = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_60s_SANO';

if ~exist(carpeta_out,'dir')
    mkdir(carpeta_out)
end

archivo_salida = fullfile( ...
    carpeta_out, ...
    'biomarcadores_sanos_60s.xlsx');

Fs_objetivo = 500;

% Ventanas dentro del tramo continuo formado por dos ventanas de 120 s.
ventanas_rel = [ ...
      0   60; ...
    120  180];

etiquetas_momento = [ ...
    "SANO_minuto_1"; ...
    "SANO_minuto_3"];

orden_momento = [1; 2];

% Controles minimos de calidad.
% En las ventanas de 60 s se exigen al menos 30 ondas P y 30 ondas T.
% Todas las ondas validas se utilizan conjuntamente para la correlacion.
min_R_60s = 20;
min_RR_validos = 20;
min_ondas = 30;

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
        R240 = R240(R240 >= 1 & R240 <= numel(ecg240));

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

            B = calcular_biomarcadores_RS_60s( ...
                ecg60, ...
                R60, ...
                Fs_objetivo, ...
                min_RR_validos, ...
                min_ondas);

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
                B.DF_residual_Hz, ...
                B.P_NumOndas, ...
                B.P_CorrIntraMedia, ...
                B.P_CorrIntraStd, ...
                B.P_AmpMedia, ...
                B.P_AmpStd, ...
                B.P_StdMedia, ...
                B.T_NumOndas, ...
                B.T_CorrIntraMedia, ...
                B.T_CorrIntraStd, ...
                B.T_AmpMedia, ...
                B.T_AmpStd, ...
                B.T_StdMedia};
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
    error('No se obtuvo ningun paciente sano con dos ventanas validas.');
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
    'DF_residual_Hz', ...
    'P_NumOndas', ...
    'P_CorrIntraMedia', ...
    'P_CorrIntraStd', ...
    'P_AmpMedia', ...
    'P_AmpStd', ...
    'P_StdMedia', ...
    'T_NumOndas', ...
    'T_CorrIntraMedia', ...
    'T_CorrIntraStd', ...
    'T_AmpMedia', ...
    'T_AmpStd', ...
    'T_StdMedia'};

T = cell2table( ...
    filas_resultados, ...
    'VariableNames',nombres_columnas);

T = sortrows(T,{'Paciente','Orden_momento'});

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

function B = calcular_biomarcadores_RS_60s( ...
    ecg,R,Fs,min_RR,min_ondas)

B = struct( ...
    'N_R',numel(R),'N_RR',0, ...
    'RR_mean',NaN,'SDNN',NaN,'RMSSD',NaN,'SDSD',NaN, ...
    'pNN20',NaN,'pNN50',NaN,'CV_RR',NaN, ...
    'SD1',NaN,'SD2',NaN,'SD1_SD2',NaN, ...
    'DF_completo_Hz',NaN,'DF_residual_Hz',NaN, ...
    'P_NumOndas',0,'P_CorrIntraMedia',NaN,'P_CorrIntraStd',NaN, ...
    'P_AmpMedia',NaN,'P_AmpStd',NaN,'P_StdMedia',NaN, ...
    'T_NumOndas',0,'T_CorrIntraMedia',NaN,'T_CorrIntraStd',NaN, ...
    'T_AmpMedia',NaN,'T_AmpStd',NaN,'T_StdMedia',NaN, ...
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

%% FILTRADO ESPECIFICO DE RS

[bRS,aRS] = butter(2,[0.5 20]/(Fs/2),'bandpass');
ecg_RS = filtfilt(bRS,aRS,ecg(:));

%% FRECUENCIA DOMINANTE

[B.DF_completo_Hz,~,~,~] = ...
    frecuencia_dominante2_RS120s(ecg_RS,Fs);

residual = cancelar_QRST_plantilla_medianaRS(ecg_RS,R,Fs);

if ~isempty(residual) && all(isfinite(residual))

    [B.DF_residual_Hz,~,~,~] = ...
        frecuencia_dominante2_RS120s(residual,Fs);
end

%% ONDA P

[ondasP,~] = extraer_ondasP_desde_R_local(ecg_RS,R,Fs);

MP = calcular_morfologia_onda_60s(ondasP,min_ondas);

B.P_NumOndas = MP.NumOndas;
B.P_CorrIntraMedia = MP.CorrIntraMedia;
B.P_CorrIntraStd = MP.CorrIntraStd;
B.P_AmpMedia = MP.AmpMedia;
B.P_AmpStd = MP.AmpStd;
B.P_StdMedia = MP.StdMedia;

%% ONDA T

[ondasT,~] = extraer_ondasT_desde_R_local(ecg_RS,R,Fs);

MT = calcular_morfologia_onda_60s(ondasT,min_ondas);

B.T_NumOndas = MT.NumOndas;
B.T_CorrIntraMedia = MT.CorrIntraMedia;
B.T_CorrIntraStd = MT.CorrIntraStd;
B.T_AmpMedia = MT.AmpMedia;
B.T_AmpStd = MT.AmpStd;
B.T_StdMedia = MT.StdMedia;

%% CONTROL FINAL

valores = [ ...
    B.RR_mean B.SDNN B.RMSSD B.SDSD ...
    B.pNN20 B.pNN50 B.CV_RR ...
    B.SD1 B.SD2 B.SD1_SD2 ...
    B.DF_completo_Hz B.DF_residual_Hz ...
    B.P_CorrIntraMedia B.P_CorrIntraStd ...
    B.P_AmpMedia B.P_AmpStd B.P_StdMedia ...
    B.T_CorrIntraMedia B.T_CorrIntraStd ...
    B.T_AmpMedia B.T_AmpStd B.T_StdMedia];

B.Valida = all(isfinite(valores));

end

function M = calcular_morfologia_onda_60s(ondas,min_ondas)

M = struct( ...
    'NumOndas',0, ...
    'CorrIntraMedia',NaN, ...
    'CorrIntraStd',NaN, ...
    'AmpMedia',NaN, ...
    'AmpStd',NaN, ...
    'StdMedia',NaN);

if isempty(ondas)
    return
end

ondas = double(ondas);
ondas = ondas(all(isfinite(ondas),2),:);

M.NumOndas = size(ondas,1);

if M.NumOndas < min_ondas
    return
end

onda_std = std(ondas,0,1,'omitnan');
amplitudes = max(ondas,[],2)-min(ondas,[],2);

M.AmpMedia = mean(amplitudes,'omitnan');
M.AmpStd = std(amplitudes,0,'omitnan');
M.StdMedia = mean(onda_std,'omitnan');

R_ondas = corrcoef(ondas');


idx_sup = triu(true(size(R_ondas)),1);
correlaciones = R_ondas(idx_sup);
correlaciones = correlaciones(isfinite(correlaciones));

if ~isempty(correlaciones)
    M.CorrIntraMedia = mean(correlaciones,'omitnan');
    M.CorrIntraStd = std(correlaciones,0,'omitnan');
end

end

function [DF, f_axis, Pxx, peak_power] = ...
    frecuencia_dominante2_RS120s(x, Fs)

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

%% Buscar pico dominante entre 0.5 y 2 Hz

idx = f_axis >= 0.5 & f_axis <= 2.0;

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

if DF <= 0.55 || DF >= 1.95
    DF = NaN;
    peak_power = NaN;
    return
end

end

function [x_residual, plantilla, t_plantilla, latidos_validos] = cancelar_QRST_plantilla_medianaRS(x, locs_R, Fs)

% CANCELAR_QRST_PLANTILLA_MEDIANARS
% Cancela los complejos QRS-T mediante sustraccion de una plantilla mediana.
%
% Version para RS:
%   - 60 ms antes del R
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

pre_R  = round(0.06 * Fs);   % 60 ms antes del R
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

function [ondasP, duraciones_ms, locs_P, inicios_P, fines_P] = extraer_ondasP_desde_R_local(x, locs_R, Fs)

% EXTRAER_ONDASP_DESDE_R_LOCAL
%
% Version robusta y permisiva para analisis morfologico de onda P.
%
% Objetivo:
%   - No perder ondas P visibles de baja amplitud.
%   - Evitar depender demasiado de findpeaks.
%   - Usar umbrales relativos al ruido local.
%   - Buscar la P en una zona fisiologica amplia antes del QRS.
%
% Uso recomendado:
%   ECG filtrado 0.5-20 Hz, sin cancelacion QRS.

x = x(:);
locs_R = limpiar_locs_local(locs_R, length(x));

ondasP_cell = cell(length(locs_R),1);
duraciones_tmp = nan(length(locs_R),1);

locs_P_tmp = nan(length(locs_R),1);
inicios_P_tmp = nan(length(locs_R),1);
fines_P_tmp = nan(length(locs_R),1);

cont = 0;

%% Longitud del segmento P guardado
preP  = round(0.07 * Fs);   % 70 ms antes del centro
postP = round(0.09 * Fs);   % 90 ms despues del centro

for i = 2:length(locs_R)

    R = locs_R(i);
    Rprev = locs_R(i-1);

    RR_prev = (R - Rprev) / Fs;

    if RR_prev < 0.30 || RR_prev > 1.5
        continue
    end

    %% ========================================================
    % 1) Zona de busqueda de P.
    % Se probo de -180 a -90 ms, pero no se detectaban todas las ondas.
    % Por eso se usa una ventana mas amplia: -300 a -50 ms.
    %% ========================================================

    t_min = max(-0.30, -0.45 * RR_prev);   % limite mas lejano
    t_max = -0.05;                          % limite mas cercano al QRS

    ini_busq = R + round(t_min * Fs);
    fin_busq = R + round(t_max * Fs);

    if ini_busq < 1 || fin_busq > length(x) || ini_busq >= fin_busq
        continue
    end

    seg = x(ini_busq:fin_busq);

    if numel(seg) < round(0.08*Fs) || any(~isfinite(seg))
        continue
    end

    %% ========================================================
    % 2) Preprocesamiento local simplificado.
    % Se centra el segmento y se aplica suavizado ligero.
    %% ========================================================

    n = numel(seg); 

    seg_dt = seg - mean(seg, 'omitnan');

    %% Suavizado ligero
    win_suave = max(3, round(0.012 * Fs));
    seg_suave = movmean(seg_dt, win_suave);

    %% ========================================================
    % 3) Estimar ruido local de forma robusta.
    %% ========================================================

    ruido = 1.4826 * median(abs(seg_suave - median(seg_suave, 'omitnan')), 'omitnan');

    if ~isfinite(ruido) || ruido <= 0
        ruido = std(seg_suave, 0, 'omitnan');
    end

    if ~isfinite(ruido) || ruido <= 0
        continue
    end

    %% ========================================================
    % 4) Buscar candidatos.
    %
    % Primero se intenta con findpeaks.
    % Si no hay picos claros, se usa maximo absoluto en la zona.
    %% ========================================================

    prom_min = max(0.004, 0.8 * ruido);
    dist_min = max(1, round(0.035 * Fs));

    candidatos_loc = [];
    candidatos_score = [];

    try
        [pks_pos, locs_pos] = findpeaks(seg_suave, ...
            'MinPeakProminence', prom_min, ...
            'MinPeakDistance', dist_min);

        if ~isempty(pks_pos)
            candidatos_loc = [candidatos_loc; locs_pos(:)];
            candidatos_score = [candidatos_score; abs(pks_pos(:))];
        end
    catch
    end

    try
        [pks_neg, locs_neg] = findpeaks(-seg_suave, ...
            'MinPeakProminence', prom_min, ...
            'MinPeakDistance', dist_min);

        if ~isempty(pks_neg)
            candidatos_loc = [candidatos_loc; locs_neg(:)];
            candidatos_score = [candidatos_score; abs(pks_neg(:))];
        end
    catch
    end

    %% Fallback: si no hay findpeaks, usar maximo absoluto
    if isempty(candidatos_loc)

        [amp_abs, idx_abs] = max(abs(seg_suave));

        if ~isfinite(amp_abs)
            continue
        end

        % Umbral permisivo pero relativo al ruido.
        if amp_abs < max(0.006, 0.9 * ruido)
            continue
        end

        locP_rel = idx_abs;

    else

        %% Preferencia: candidato con mas amplitud, evitando extremos
        margen = round(0.02 * Fs);

        idx_valid = candidatos_loc > margen & candidatos_loc < (n - margen);

        if any(idx_valid)
            candidatos_loc = candidatos_loc(idx_valid);
            candidatos_score = candidatos_score(idx_valid);
        end

        [~, idx_best] = max(candidatos_score);
        locP_rel = candidatos_loc(idx_best);
    end

    locP = ini_busq + locP_rel - 1;

    %% ========================================================
    % 5) Extraer segmento P.
    %% ========================================================

    iniP = locP - preP;
    finP = locP + postP;

    % Evitar que entre QRS.
    finP_max = R - round(0.025 * Fs);

    if finP > finP_max
        finP = finP_max;
    end

    if iniP < 1 || finP > length(x) || iniP >= finP
        continue
    end

    p = x(iniP:finP);

    if any(~isfinite(p))
        continue
    end

    %% ========================================================
    % 6) Normalizar longitud del segmento.
    %% ========================================================

    L_obj = preP + postP + 1;

    if length(p) ~= L_obj
        t_old = linspace(0, 1, length(p));
        t_new = linspace(0, 1, L_obj);
        p = interp1(t_old, p, t_new, 'linear', 'extrap')';
    end

    %% Baseline local
    n_base = max(1, round(0.02 * Fs));
    baseline = mean(p(1:min(n_base, length(p))), 'omitnan');
    p = p - baseline;

    %% ========================================================
    % 7) Criterios de calidad suaves.
    % No se exige una P perfecta; solo se descartan segmentos planos
    % o claramente anormales.
    %% ========================================================

    ampP = max(p, [], 'omitnan') - min(p, [], 'omitnan');
    stdP = std(p, 0, 'omitnan');

    if ~isfinite(ampP) || ~isfinite(stdP)
        continue
    end

    if stdP < 0.002
        continue
    end

    if ampP < 0.008
        continue
    end

    %% Evitar segmentos dominados por salto brusco tipo QRS
    dp = abs(diff(p));
    
    if max(dp, [], 'omitnan') > 0.30
        continue
    end

    %% ========================================================
    % 8) Inicio y fin aproximados de P.
    %% ========================================================

    y = p - median(p, 'omitnan');
    yabs = abs(y);

    pico = max(yabs, [], 'omitnan');

    if ~isfinite(pico) || pico <= 0
        continue
    end

    umbral = 0.12 * pico;
    idx_sup = find(yabs >= umbral);

    if isempty(idx_sup)
        continue
    end

    ini_rel = idx_sup(1);
    fin_rel = idx_sup(end);

    dur_ms = 1000 * (fin_rel - ini_rel + 1) / Fs;

    %% Criterio amplio
    if dur_ms < 15 || dur_ms > 240
        continue
    end

    %% Guardar

    cont = cont + 1;

    ondasP_cell{cont} = p(:)';
    duraciones_tmp(cont) = dur_ms;

    locs_P_tmp(cont) = locP;
    inicios_P_tmp(cont) = iniP + ini_rel - 1;
    fines_P_tmp(cont) = iniP + fin_rel - 1;

end

%% Salidas

if cont == 0

    ondasP = [];
    duraciones_ms = [];
    locs_P = [];
    inicios_P = [];
    fines_P = [];

else

    ondasP = cell2mat(ondasP_cell(1:cont));
    duraciones_ms = duraciones_tmp(1:cont);

    locs_P = locs_P_tmp(1:cont);
    inicios_P = inicios_P_tmp(1:cont);
    fines_P = fines_P_tmp(1:cont);

end

end

function [ondasT, duraciones_ms, locs_T, inicios_T, fines_T] = extraer_ondasT_desde_R_local(x, locs_R, Fs)

% EXTRAER_ONDAST_DESDE_R_LOCAL
%
% Version robusta para analisis morfologico de onda T.
%
% Objetivo:
%   - Detectar ondas T anchas o poco picudas.
%   - No depender exclusivamente de findpeaks.
%   - Usar una ventana fisiologica adaptada al RR.
%   - Evitar contaminarse con el siguiente QRS.
%
% Uso recomendado:
%   ECG filtrado 0.5-20 Hz, sin cancelacion QRS.

x = x(:);
locs_R = limpiar_locs_local(locs_R, length(x));

ondasT_cell = cell(length(locs_R), 1);
duraciones_tmp = nan(length(locs_R), 1);

locs_T_tmp = nan(length(locs_R), 1);
inicios_T_tmp = nan(length(locs_R), 1);
fines_T_tmp = nan(length(locs_R), 1);

cont = 0;

%% Longitud del segmento T guardado
preT  = round(0.10 * Fs);   % 100 ms antes del centro
postT = round(0.14 * Fs);   % 140 ms despues del centro

for i = 1:length(locs_R)-1

    R = locs_R(i);
    Rnext = locs_R(i+1);

    RR_next = (Rnext - R) / Fs;

    if RR_next < 0.30 || RR_next > 1.5
        continue
    end

    %% ========================================================
    % 1) Zona de busqueda de T.
    %
    % Empieza despues del QRS y termina antes del siguiente QRS.
    % Se adapta al RR para no invadir el siguiente latido.
    %% ========================================================

    t_min = 0.10;
    t_max = min(0.65 * RR_next, 0.48);

    ini_busq = R + round(t_min * Fs);
    fin_busq = R + round(t_max * Fs);

    % No acercarse demasiado al siguiente QRS.
    fin_busq = min(fin_busq, Rnext - round(0.04 * Fs));

    if ini_busq < 1 || fin_busq > length(x) || ini_busq >= fin_busq
        continue
    end

    seg = x(ini_busq:fin_busq);

    if numel(seg) < round(0.10 * Fs) || any(~isfinite(seg))
        continue
    end

    %% ========================================================
    % 2) Preprocesamiento local simplificado.
    % Se centra el segmento y se aplica suavizado ligero.
    %% ========================================================

    n = numel(seg); 

    seg_dt = seg - mean(seg, 'omitnan');

    %% Suavizado algo mayor que en P porque la T es mas ancha
    win_suave = max(3, round(0.025 * Fs));
    seg_suave = movmean(seg_dt, win_suave);

    %% ========================================================
    % 3) Ruido local robusto.
    %% ========================================================

    ruido = 1.4826 * median(abs(seg_suave - median(seg_suave, 'omitnan')), 'omitnan');

    if ~isfinite(ruido) || ruido <= 0
        ruido = std(seg_suave, 0, 'omitnan');
    end

    if ~isfinite(ruido) || ruido <= 0
        continue
    end

    %% ========================================================
    % 4) Buscar candidatos positivos y negativos.
    %
    % Si findpeaks falla, se usa maximo absoluto.
    %% ========================================================

    prom_min = max(0.006, 0.7 * ruido);
    dist_min = max(1, round(0.06 * Fs));

    candidatos_loc = [];
    candidatos_score = [];

    try
        [pks_pos, locs_pos] = findpeaks(seg_suave, ...
            'MinPeakProminence', prom_min, ...
            'MinPeakDistance', dist_min);

        if ~isempty(pks_pos)
            candidatos_loc = [candidatos_loc; locs_pos(:)];
            candidatos_score = [candidatos_score; abs(pks_pos(:))];
        end
    catch
    end

    try
        [pks_neg, locs_neg] = findpeaks(-seg_suave, ...
            'MinPeakProminence', prom_min, ...
            'MinPeakDistance', dist_min);

        if ~isempty(pks_neg)
            candidatos_loc = [candidatos_loc; locs_neg(:)];
            candidatos_score = [candidatos_score; abs(pks_neg(:))];
        end
    catch
    end

    %% Fallback para T ancha y suave
    if isempty(candidatos_loc)

        [amp_abs, idx_abs] = max(abs(seg_suave));

        if ~isfinite(amp_abs)
            continue
        end

        if amp_abs < max(0.008, 0.8 * ruido)
            continue
        end

        locT_rel = idx_abs;

    else

        %% Evitar candidatos pegados a los bordes de la ventana
        margen = round(0.025 * Fs);

        idx_valid = candidatos_loc > margen & candidatos_loc < (n - margen);

        if any(idx_valid)
            candidatos_loc = candidatos_loc(idx_valid);
            candidatos_score = candidatos_score(idx_valid);
        end

        [~, idx_best] = max(candidatos_score);
        locT_rel = candidatos_loc(idx_best);
    end

    locT = ini_busq + locT_rel - 1;

    %% ========================================================
    % 5) Extraer segmento T.
    %% ========================================================

    iniT = locT - preT;
    finT = locT + postT;

    % Evitar el siguiente QRS.
    finT_max = Rnext - round(0.04 * Fs);

    if finT > finT_max
        finT = finT_max;
    end

    if iniT < 1 || finT > length(x) || iniT >= finT
        continue
    end

    t = x(iniT:finT);

    if any(~isfinite(t))
        continue
    end

    %% ========================================================
    % 6) Normalizar longitud del segmento.
    %% ========================================================

    L_obj = preT + postT + 1;

    if length(t) ~= L_obj
        t_old = linspace(0, 1, length(t));
        t_new = linspace(0, 1, L_obj);
        t = interp1(t_old, t, t_new, 'linear', 'extrap')';
    end

    %% Baseline local
    n_base = max(1, round(0.03 * Fs));
    baseline = mean(t(1:min(n_base, length(t))), 'omitnan');
    t = t - baseline;

    %% ========================================================
    % 7) Criterios de calidad suaves.
    %% ========================================================

    ampT = max(t, [], 'omitnan') - min(t, [], 'omitnan');
    stdT = std(t, 0, 'omitnan');

    if ~isfinite(ampT) || ~isfinite(stdT)
        continue
    end

    if stdT < 0.003
        continue
    end

    if ampT < 0.012
        continue
    end

    %% Evitar segmentos dominados por QRS residual o saltos muy bruscos
    dt = abs(diff(t));

    if max(dt, [], 'omitnan') > 0.40
        continue
    end
   

    %% ========================================================
    % 8) Inicio y fin aproximados de T.
    %% ========================================================

    y = t - median(t, 'omitnan');
    yabs = abs(y);

    pico = max(yabs, [], 'omitnan');

    if ~isfinite(pico) || pico <= 0
        continue
    end

    umbral = 0.10 * pico;
    idx_sup = find(yabs >= umbral);

    if isempty(idx_sup)
        continue
    end

    ini_rel = idx_sup(1);
    fin_rel = idx_sup(end);

    dur_ms = 1000 * (fin_rel - ini_rel + 1) / Fs;

    if dur_ms < 40 || dur_ms > 360
        continue
    end

    %% Guardar

    cont = cont + 1;

    ondasT_cell{cont} = t(:)';
    duraciones_tmp(cont) = dur_ms;

    locs_T_tmp(cont) = locT;
    inicios_T_tmp(cont) = iniT + ini_rel - 1;
    fines_T_tmp(cont) = iniT + fin_rel - 1;

end

%% Salidas

if cont == 0

    ondasT = [];
    duraciones_ms = [];
    locs_T = [];
    inicios_T = [];
    fines_T = [];

else

    ondasT = cell2mat(ondasT_cell(1:cont));
    duraciones_ms = duraciones_tmp(1:cont);

    locs_T = locs_T_tmp(1:cont);
    inicios_T = inicios_T_tmp(1:cont);
    fines_T = fines_T_tmp(1:cont);

end

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