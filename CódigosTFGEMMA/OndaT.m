clear
clc
close all

%% ============================================================
% CÁLCULO DE BIOMARCADORES MORFOLÓGICOS DE LA ONDA T
%
% Este script calcula biomarcadores de la onda T en ventanas de ritmo
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
%   - Extrae ondas T con la función local extraer_ondasT_desde_R_local.
%   - Calcula:
%       NumOndasT
%       CorrIntraMedia
%       CorrIntraStd
%       AmpMedia
%       AmpStd
%       StdMedia
%
% NOTA:
%   - Se eliminan las duraciones de onda T como biomarcador principal,
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

carpeta_out = 'C:\Users\Emma\Documents\MATLAB\analisis_ondaTHIBRIDOnuevo';

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

        ondasT = [];

        n_ventanas_totales = numel(archivos_reg);
        n_ventanas_usadas = 0;
        n_descartadas_sin_R = 0;
        n_descartadas_pocas_R = 0;
        n_descartadas_senal = 0;
        n_ondasT_total = 0;

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

            %% Filtrado para onda T

            [bT, aT] = butter(2, [0.5 20] / (Fs / 2), 'bandpass');
            xT = filtfilt(bT, aT, x);

            if isempty(xT) || any(~isfinite(xT))
                n_descartadas_senal = n_descartadas_senal + 1;
                continue
            end

            %% Limpiar locs_R guardados
            % No se vuelven a ajustar. Ya vienen ajustados del dataset híbrido.

            locs_R = unique(round(locs_R(:)));
            locs_R = locs_R(isfinite(locs_R));
            locs_R = locs_R(locs_R >= 1 & locs_R <= length(xT));

            if isempty(locs_R) || numel(locs_R) < 3
                n_descartadas_pocas_R = n_descartadas_pocas_R + 1;
                continue
            end

            %% Extraer ondas T

            [ondasT_vent, ~] = extraer_ondasT_desde_R_local(xT, locs_R, Fs);

            if ~isempty(ondasT_vent)

                ondasT = [ondasT; ondasT_vent]; %#ok<AGROW>

                n_ventanas_usadas = n_ventanas_usadas + 1;
                n_ondasT_total = n_ondasT_total + size(ondasT_vent, 1);
            end
        end

        %% ====================================================
        % CONTROL MÍNIMO DE ONDAS T
        %% ====================================================

        if isempty(ondasT)

            fprintf('Sin ondas T válidas en %s\n', reg);
            fprintf('  Ventanas totales       : %d\n', n_ventanas_totales);
            fprintf('  Ventanas usadas        : %d\n', n_ventanas_usadas);
            fprintf('  Sin locs_R             : %d\n', n_descartadas_sin_R);
            fprintf('  Pocas R                : %d\n', n_descartadas_pocas_R);
            fprintf('  Señal no válida        : %d\n', n_descartadas_senal);

            continue
        end

        if size(ondasT, 1) < 30

            fprintf('Pocas ondas T en %s: %d\n', reg, size(ondasT,1));
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

        % Onda T media y variabilidad punto a punto
        T_media = mean(ondasT, 1, 'omitnan');
        T_std   = std(ondasT, 0, 1, 'omitnan');

        num_ondasT = size(ondasT, 1);

        %% Amplitud pico a pico

        amplitudes = max(ondasT, [], 2) - min(ondasT, [], 2);

        amp_media = mean(amplitudes, 'omitnan');
        amp_std   = std(amplitudes, 0, 'omitnan');

        %% Variabilidad morfológica media

        std_media = mean(T_std, 'omitnan');

        %% ====================================================
        % CORRELACIÓN INTRAPACIENTE POR BLOQUES DE 10 ONDAS T
        %% ====================================================

        n_bloques = floor(size(ondasT, 1) / 10);
        correlaciones_paciente = [];

        for b = 1:n_bloques

            bloque = ondasT((b - 1) * 10 + 1 : b * 10, :);

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
            num_ondasT, ...
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
            'NumOndasT', ...
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

        resultados_grupo(idx_res).T_media = T_media;
        resultados_grupo(idx_res).T_std = T_std;

        resultados_grupo(idx_res).NumOndasT = num_ondasT;

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
        fprintf('  Ondas T válidas       : %d\n', num_ondasT);
        fprintf('  Ventanas totales      : %d\n', n_ventanas_totales);
        fprintf('  Ventanas usadas       : %d\n', n_ventanas_usadas);
        fprintf('  Corr intra media      : %.3f\n', corr_intra_media);
        fprintf('  Amp T media           : %.6f\n', amp_media);
        fprintf('  Std T media           : %.6f\n', std_media);

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

archivo_excel = fullfile(carpeta_out, 'features_ondaT.xlsx');
archivo_mat   = fullfile(carpeta_out, 'ondasT_features_completas.mat');

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

fprintf('\nFIN ANALISIS ONDA T - HIBRIDO_R\n');


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

function [ondasT, duraciones_ms, locs_T, inicios_T, fines_T] = extraer_ondasT_desde_R_local(x, locs_R, Fs)

% EXTRAER_ONDAST_DESDE_R_LOCAL
%
% Versión robusta para análisis morfológico de onda T.
%
% Objetivo:
%   - Detectar ondas T anchas o poco picudas.
%   - No depender exclusivamente de findpeaks.
%   - Usar una ventana fisiológica adaptada al RR.
%   - Evitar contaminarse con el siguiente QRS.
%
% Uso recomendado:
%   ECG filtrado 0.5-20 Hz, sin cancelación QRS.

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
postT = round(0.14 * Fs);   % 140 ms después del centro

for i = 1:length(locs_R)-1

    R = locs_R(i);
    Rnext = locs_R(i+1);

    RR_next = (Rnext - R) / Fs;

    if RR_next < 0.30 || RR_next > 1.5
        continue
    end

    %% ========================================================
    % 1) Zona de búsqueda de T
    %
    % Empieza después del QRS y termina antes del siguiente QRS.
    % Se adapta al RR para no invadir el siguiente latido.
    %% ========================================================

    t_min = 0.10;
    t_max = min(0.65 * RR_next, 0.48);

    ini_busq = R + round(t_min * Fs);
    fin_busq = R + round(t_max * Fs);

    % No acercarse demasiado al siguiente QRS
    fin_busq = min(fin_busq, Rnext - round(0.04 * Fs));

    if ini_busq < 1 || fin_busq > length(x) || ini_busq >= fin_busq
        continue
    end

    seg = x(ini_busq:fin_busq);

    if numel(seg) < round(0.10 * Fs) || any(~isfinite(seg))
        continue
    end

    %% ========================================================
    % 2) Preprocesamiento local simplificado
    % Se centra el segmento y se aplica suavizado ligero.
    %% ========================================================

    n = numel(seg); 

    seg_dt = seg - mean(seg, 'omitnan');

    %% Suavizado algo mayor que en P porque la T es más ancha
    win_suave = max(3, round(0.025 * Fs));
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
    % 4) Buscar candidatos positivos y negativos
    %
    % Si findpeaks falla, se usa el máximo absoluto de la zona.
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
    % 5) Extraer segmento T
    %% ========================================================

    iniT = locT - preT;
    finT = locT + postT;

    % Evitar el siguiente QRS
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
    % 6) Normalizar longitud del segmento
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
    % 7) Criterios de calidad suaves
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
    % 8) Inicio y fin aproximados de T
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