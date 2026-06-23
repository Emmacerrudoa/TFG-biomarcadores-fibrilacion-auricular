clear
clc
close all

%% ============================================================
% CÁLCULO DE BIOMARCADORES MORFOLÓGICOS DE LA ONDA P
%
% Este script calcula biomarcadores de la onda P en ventanas de ritmo
% sinusal procedentes del dataset generado previamente.
%
% Dataset:
%   dataset_final_HIBRIDO
%
% IMPORTANTE:
%   - Los locs_R ya están guardados y ajustados en el preprocesado híbrido.
%   - Este script NO vuelve a detectar R.
%   - Este script NO vuelve a aplicar ajustar_R_al_pico.
%
% Grupos analizados:
%   - SANO
%   - FA_PAROXISTICA_RS
%
% Para cada paciente:
%   - Extrae ondas P con la función local extraer_ondasP_desde_R_local.
%   - Calcula:
%       NumOndasP
%       CorrIntraMedia
%       CorrIntraStd
%       AmpMedia
%       AmpStd
%       StdMedia
%
% NOTA:
%   - Se eliminan las duraciones de onda P como biomarcador principal,
%     porque dependen de una delimitación por umbral dentro de una ventana fija.
%
% Requisitos:
%   - Archivos .mat generados en la fase de segmentación.
%   - MATLAB con Signal Processing Toolbox.
%% ============================================================

%% CONFIGURACIÓN GENERAL

% Modificar estas rutas según la ubicación local del dataset generado.
grupos_proc = { ...
    struct('ruta','C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\SANO', ...
           'grupo','SANO'), ...
    struct('ruta','C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\FA_PAROXISTICA_RS', ...
           'grupo','FA_PAROXISTICA_RS') ...
};

carpeta_out = 'C:\Users\Emma\Documents\MATLAB\analisis_ondaPHIBRIDOnuevo';

if ~exist(carpeta_out, 'dir')
    mkdir(carpeta_out);
end

Fs_default = 500;

tabla_features = table;

resultados_sanos = struct([]);
resultados_parox = struct([]);

correlaciones_intra_sanos = [];
correlaciones_intra_parox = [];

%% ============================================================
% PROCESAR GRUPOS
%% ============================================================

for gg = 1:numel(grupos_proc)

    carpeta_grupo = grupos_proc{gg}.ruta;
    nombre_grupo  = grupos_proc{gg}.grupo;

    fprintf('\n============================================\n');
    fprintf('Procesando grupo: %s\n', nombre_grupo);
    fprintf('Carpeta: %s\n', carpeta_grupo);
    fprintf('============================================\n');

    archivos = dir(fullfile(carpeta_grupo, '*.mat'));

    if isempty(archivos)
        warning('No hay archivos .mat en %s', carpeta_grupo);
        continue
    end

    %% Obtener nombres de registro

    registros = cell(numel(archivos), 1);

    for k = 1:numel(archivos)

        ruta_archivo = fullfile(carpeta_grupo, archivos(k).name);
        D = load(ruta_archivo, 'nombre_registro');

        if isfield(D, 'nombre_registro') && ~isempty(D.nombre_registro)
            registros{k} = char(string(D.nombre_registro));
        else
            registros{k} = '';
        end
    end

    registros_unicos = unique(registros);
    registros_unicos(cellfun(@isempty, registros_unicos)) = [];

    %% ========================================================
    % PROCESAR PACIENTES / REGISTROS
    %% ========================================================

    resultados_grupo = struct([]);
    correlaciones_intra_grupo = [];

    for r = 1:numel(registros_unicos)

        reg = registros_unicos{r};

        fprintf('\n--------------------------------------------\n');
        fprintf('Procesando %s - %s\n', nombre_grupo, reg);
        fprintf('--------------------------------------------\n');

        idx_reg = strcmp(registros, reg);
        archivos_reg = archivos(idx_reg);

        ondasP = [];

        n_ventanas_totales = numel(archivos_reg);
        n_ventanas_usadas = 0;
        n_descartadas_sin_R = 0;
        n_descartadas_pocas_R = 0;
        n_descartadas_senal = 0;
        n_ondasP_total = 0;

        %% ====================================================
        % PROCESAR VENTANAS DEL PACIENTE
        %% ====================================================

        for k = 1:numel(archivos_reg)

            ruta_archivo = fullfile(carpeta_grupo, archivos_reg(k).name);

            D = load(ruta_archivo, ...
                'ventana', 'locs_R', 'Fs', ...
                'nombre_registro', 'ID_global', 't_ini', 't_fin', ...
                'ritmo_ventana', 'tipo_registro');

            if ~isfield(D, 'ventana') || isempty(D.ventana)
                n_descartadas_senal = n_descartadas_senal + 1;
                continue
            end

            if ~isfield(D, 'locs_R') || isempty(D.locs_R)
                n_descartadas_sin_R = n_descartadas_sin_R + 1;
                continue
            end

            x = double(D.ventana(:));
            locs_R = round(D.locs_R(:));

            if isfield(D, 'Fs') && ~isempty(D.Fs) && isfinite(D.Fs)
                Fs = D.Fs;
            else
                Fs = Fs_default;
            end

            if isempty(x) || any(~isfinite(x))
                n_descartadas_senal = n_descartadas_senal + 1;
                continue
            end

            %% Filtrado para onda P

            [bP, aP] = butter(2, [0.5 20] / (Fs / 2), 'bandpass');
            xP = filtfilt(bP, aP, x);

            if isempty(xP) || any(~isfinite(xP))
                n_descartadas_senal = n_descartadas_senal + 1;
                continue
            end

            %% Limpiar locs_R guardados
            % No se vuelven a ajustar. Ya vienen ajustados del dataset híbrido.

            locs_R = limpiar_locs_local(locs_R, length(xP));

            if isempty(locs_R) || numel(locs_R) < 3
                n_descartadas_pocas_R = n_descartadas_pocas_R + 1;
                continue
            end

            %% Extraer ondas P
            % La función devuelve también duraciones, pero aquí se ignoran.

            [ondasP_vent, ~] = extraer_ondasP_desde_R_local(xP, locs_R, Fs);

            if ~isempty(ondasP_vent)

                ondasP = [ondasP; ondasP_vent]; %#ok<AGROW>

                n_ventanas_usadas = n_ventanas_usadas + 1;
                n_ondasP_total = n_ondasP_total + size(ondasP_vent, 1);
            end
        end

        %% ====================================================
        % CONTROL MÍNIMO DE ONDAS P
        %% ====================================================

        if isempty(ondasP)

            fprintf('Sin ondas P válidas en %s\n', reg);
            fprintf('  Ventanas totales       : %d\n', n_ventanas_totales);
            fprintf('  Ventanas usadas        : %d\n', n_ventanas_usadas);
            fprintf('  Sin locs_R             : %d\n', n_descartadas_sin_R);
            fprintf('  Pocas R                : %d\n', n_descartadas_pocas_R);
            fprintf('  Señal no válida        : %d\n', n_descartadas_senal);

            continue
        end

        if size(ondasP, 1) < 30

            fprintf('Pocas ondas P en %s: %d\n', reg, size(ondasP,1));
            fprintf('  Ventanas totales       : %d\n', n_ventanas_totales);
            fprintf('  Ventanas usadas        : %d\n', n_ventanas_usadas);
            fprintf('  Sin locs_R             : %d\n', n_descartadas_sin_R);
            fprintf('  Pocas R                : %d\n', n_descartadas_pocas_R);
            fprintf('  Señal no válida        : %d\n', n_descartadas_senal);

            continue
        end

        %% ====================================================
        % BIOMARCADORES DEL PACIENTE
        %% ====================================================

        % Onda P media y variabilidad punto a punto
        P_media = mean(ondasP, 1, 'omitnan');
        P_std   = std(ondasP, 0, 1, 'omitnan');

        num_ondasP = size(ondasP, 1);

        %% Amplitud pico a pico

        amplitudes = max(ondasP, [], 2) - min(ondasP, [], 2);

        amp_media = mean(amplitudes, 'omitnan');
        amp_std   = std(amplitudes, 0, 'omitnan');

        %% Variabilidad morfológica media

        std_media = mean(P_std, 'omitnan');

        %% ====================================================
        % CORRELACIÓN INTRAPACIENTE POR BLOQUES DE 10 ONDAS P
        %% ====================================================

        n_bloques = floor(size(ondasP, 1) / 10);
        correlaciones_paciente = [];

        for b = 1:n_bloques

            bloque = ondasP((b - 1) * 10 + 1 : b * 10, :);

            if size(bloque, 1) < 2
                continue
            end

            R_bloque = corrcoef(bloque');

            idx_sup = triu(true(size(R_bloque)), 1);
            correlacion_bloque = R_bloque(idx_sup);

            correlacion_bloque = correlacion_bloque(isfinite(correlacion_bloque));

            correlaciones_paciente = [correlaciones_paciente; correlacion_bloque(:)]; %#ok<AGROW>
        end

        if isempty(correlaciones_paciente)
            corr_intra_media = NaN;
            corr_intra_std   = NaN;
        else
            corr_intra_media = mean(correlaciones_paciente, 'omitnan');
            corr_intra_std   = std(correlaciones_paciente, 0, 'omitnan');
        end

        %% ====================================================
        % AÑADIR FILA A LA TABLA
        %% ====================================================

        nueva_fila = table( ...
            {reg}, ...
            {nombre_grupo}, ...
            num_ondasP, ...
            corr_intra_media, ...
            corr_intra_std, ...
            amp_media, ...
            amp_std, ...
            std_media, ...
            n_ventanas_totales, ...
            n_ventanas_usadas, ...
            n_descartadas_sin_R, ...
            n_descartadas_pocas_R, ...
            n_descartadas_senal, ...
            'VariableNames', { ...
            'Paciente', ...
            'Grupo', ...
            'NumOndasP', ...
            'CorrIntraMedia', ...
            'CorrIntraStd', ...
            'AmpMedia', ...
            'AmpStd', ...
            'StdMedia', ...
            'N_ventanas_totales', ...
            'N_ventanas_usadas', ...
            'N_descartadas_sin_R', ...
            'N_descartadas_pocas_R', ...
            'N_descartadas_senal'} ...
        );

        tabla_features = [tabla_features; nueva_fila]; %#ok<AGROW>

        correlaciones_intra_grupo = [correlaciones_intra_grupo; correlaciones_paciente(:)]; %#ok<AGROW>

        %% ====================================================
        % GUARDAR RESULTADOS POR PACIENTE
        %% ====================================================

        idx_res = numel(resultados_grupo) + 1;

        resultados_grupo(idx_res).registro = reg;
        resultados_grupo(idx_res).grupo = nombre_grupo;

        resultados_grupo(idx_res).P_media = P_media;
        resultados_grupo(idx_res).P_std = P_std;

        resultados_grupo(idx_res).NumOndasP = num_ondasP;

        resultados_grupo(idx_res).CorrIntraMedia = corr_intra_media;
        resultados_grupo(idx_res).CorrIntraStd = corr_intra_std;

        resultados_grupo(idx_res).AmpMedia = amp_media;
        resultados_grupo(idx_res).AmpStd = amp_std;

        resultados_grupo(idx_res).StdMedia = std_media;

        resultados_grupo(idx_res).N_ventanas_totales = n_ventanas_totales;
        resultados_grupo(idx_res).N_ventanas_usadas = n_ventanas_usadas;
        resultados_grupo(idx_res).N_descartadas_sin_R = n_descartadas_sin_R;
        resultados_grupo(idx_res).N_descartadas_pocas_R = n_descartadas_pocas_R;
        resultados_grupo(idx_res).N_descartadas_senal = n_descartadas_senal;

        %% Mostrar resumen por paciente

        fprintf('Paciente %s procesado correctamente\n', reg);
        fprintf('  Ondas P válidas       : %d\n', num_ondasP);
        fprintf('  Ventanas totales      : %d\n', n_ventanas_totales);
        fprintf('  Ventanas usadas       : %d\n', n_ventanas_usadas);
        fprintf('  Corr intra media      : %.3f\n', corr_intra_media);
        fprintf('  Amp P media           : %.6f\n', amp_media);
        fprintf('  Std P media           : %.6f\n', std_media);

    end

    %% Guardar resultados del grupo en variables finales

    if strcmp(nombre_grupo, 'SANO')
        resultados_sanos = resultados_grupo;
        correlaciones_intra_sanos = correlaciones_intra_grupo;
    elseif strcmp(nombre_grupo, 'FA_PAROXISTICA_RS')
        resultados_parox = resultados_grupo;
        correlaciones_intra_parox = correlaciones_intra_grupo;
    end
end

%% ============================================================
% GUARDAR RESULTADOS
%% ============================================================

archivo_excel = fullfile(carpeta_out, 'features_ondaP.xlsx');
archivo_mat   = fullfile(carpeta_out, 'ondasP_features_completas.mat');

writetable(tabla_features, archivo_excel);

save(archivo_mat, ...
    'tabla_features', ...
    'resultados_sanos', ...
    'resultados_parox', ...
    'correlaciones_intra_sanos', ...
    'correlaciones_intra_parox', ...
    '-v7.3');

fprintf('\nTabla de features guardada en:\n%s\n', archivo_excel);
fprintf('Datos completos guardados en:\n%s\n', archivo_mat);

fprintf('\nFIN ANALISIS ONDA P - HIBRIDO_R\n');

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

function [ondasP, duraciones_ms, locs_P, inicios_P, fines_P] = extraer_ondasP_desde_R_local(x, locs_R, Fs)

% EXTRAER_ONDASP_DESDE_R_LOCAL
%
% Versión robusta/permisiva para análisis morfológico de onda P.
%
% Objetivo:
%   - No perder ondas P visibles de baja amplitud.
%   - Evitar depender demasiado de findpeaks.
%   - Usar umbrales relativos al ruido local.
%   - Buscar la P en una zona fisiológica amplia antes del QRS.
%
% Uso recomendado:
%   ECG filtrado 0.5-20 Hz, sin cancelación QRS.

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
postP = round(0.09 * Fs);   % 90 ms después del centro

for i = 2:length(locs_R)

    R = locs_R(i);
    Rprev = locs_R(i-1);

    RR_prev = (R - Rprev) / Fs;

    if RR_prev < 0.30 || RR_prev > 1.5
        continue
    end

    %% ========================================================
    % 1) Zona de búsqueda de P
    %
    % Se probó una ventana de -180 a -90 ms, pero no se detectaban todas
    % las ondas P. Por ello, se utiliza una ventana más amplia: -300 a -50 ms.
    %% ========================================================

    t_min = max(-0.30, -0.45 * RR_prev);   % límite más lejano
    t_max = -0.05;                          % límite más cercano al QRS

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
    % 2) Preprocesamiento local simplificado
    % Se centra el segmento y se aplica suavizado ligero.
    %% ========================================================

    n = numel(seg); 

    seg_dt = seg - mean(seg, 'omitnan');

    %% Suavizado ligero
    win_suave = max(3, round(0.012 * Fs));
    seg_suave = movmean(seg_dt, win_suave);

    %% ========================================================
    % 3) Estimar ruido local de forma robusta
    %% ========================================================

    ruido = 1.4826 * median(abs(seg_suave - median(seg_suave, 'omitnan')), 'omitnan');

    if ~isfinite(ruido) || ruido <= 0
        ruido = std(seg_suave, 0, 'omitnan');
    end

    if ~isfinite(ruido) || ruido <= 0
        continue
    end

    %% ========================================================
    % 4) Buscar candidatos
    %
    % Primero se intenta con findpeaks.
    % Si no hay picos claros, se usa máximo absoluto en la zona.
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

    %% Fallback: si no hay findpeaks, usar máximo absoluto
    if isempty(candidatos_loc)

        [amp_abs, idx_abs] = max(abs(seg_suave));

        if ~isfinite(amp_abs)
            continue
        end

        % Umbral muy permisivo pero relativo al ruido
        if amp_abs < max(0.006, 0.9 * ruido)
            continue
        end

        locP_rel = idx_abs;

    else

        %% Preferencia: candidato con más amplitud, pero evitando extremos
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
    % 5) Extraer segmento P
    %% ========================================================

    iniP = locP - preP;
    finP = locP + postP;

    % Evitar que entre QRS
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
    % 6) Normalizar longitud del segmento
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
    % 7) Criterios de calidad muy suaves
    % Aquí no se exige una P perfecta. Solo se descartan segmentos planos o absurdos.
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
    % 8) Inicio y fin aproximados de P
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