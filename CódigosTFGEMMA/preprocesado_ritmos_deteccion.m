clear
clc
close all

%% ============================================================
% GENERACIÓN DEL DATASET FINAL MEDIANTE DETECCIÓN HÍBRIDA DE PICOS R
%
% Este script recorre las bases de datos utilizadas en el TFG, lee la señal
% ECG y sus anotaciones mediante WFDB, remuestrea la señal a 500 Hz, detecta
% segmentos con artefactos y extrae ventanas limpias de 120 s.
%
% La localización final de los picos R se obtiene mediante un método híbrido:
%
%   - En sujetos sanos se utilizan las anotaciones beat-to-beat como
%     referencia principal.
%
%   - En registros con FA se utilizan las anotaciones .qrs como referencia
%     principal.
%
%   - La función DetectarPICOSR.m revisa y ajusta las posiciones de los
%     picos R, añadiendo candidatos solo cuando permiten completar huecos
%     largos entre latidos.
%
% Se guardan archivos .mat con:
%   ventana, locs_R, Fs, nombre_registro, nombre_base, s, w, ID_global,
%   t_ini, t_fin, tipo_registro y ritmo_ventana.
%
% Requisitos:
%   - MATLAB con Signal Processing Toolbox.
%   - Python configurado desde MATLAB.
%   - Paquetes Python: wfdb y numpy.
%   - DetectarPICOSR.m en la misma carpeta o en el path de MATLAB.
%% ============================================================

%% 0) CONFIGURACIÓN GENERAL

% Rutas de las bases de datos. Modificar estas rutas según la ubicación
% local de los archivos en cada equipo.
bases = { ...
    struct('carpeta','C:\Users\Emma\Documents\MATLAB\long-term-af-database-1.0.0\files', 'nombre','BASE_1'), ...
    struct('carpeta','C:\Users\Emma\Documents\MATLAB\mit-bih-atrial-fibrillation-database-1.0.0\files', 'nombre','BASE_2'), ...
    struct('carpeta','C:\Users\Emma\Documents\MATLAB\mit-bih-normal-sinus-rhythm-database-1.0.0', 'nombre','BASE_3_SANOS') ...
    };

% Carpeta donde se guardarán las ventanas finales y los archivos resumen.
carpeta_salida = 'C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO';

subcarpetas = {'SANO','FA_PERSISTENTE','FA_PAROXISTICA_RS','FA_PAROXISTICA_FA'};

if ~exist(carpeta_salida,'dir')
    mkdir(carpeta_salida);
end

for kk = 1:numel(subcarpetas)
    ruta_sub = fullfile(carpeta_salida, subcarpetas{kk});
    if ~exist(ruta_sub,'dir')
        mkdir(ruta_sub);
    end
end

%% 1) COMPROBAR PYTHON

try
    py.importlib.import_module('wfdb');
    py.importlib.import_module('numpy');
catch ME
    error(['MATLAB no puede importar wfdb/numpy desde Python. Mensaje: ' ME.message]);
end

%% 2) PARÁMETROS

Fs_new = 500;
dur_artefacto = 2;
dur_vent_final = 120;

resumen_global = {};
errores_global = {};
descartadas_R_global = {};

ID_global = buscar_ultimo_id_global(carpeta_salida);

fprintf('Inicio del procesamiento global\n');
fprintf('ID inicial: %d\n', ID_global);

%% 3) RECORRER BASES

for bb = 1:numel(bases)

    carpeta_base = bases{bb}.carpeta;
    nombre_base = bases{bb}.nombre;

    fprintf('\n============================================\n');
    fprintf('Procesando base: %s\n', nombre_base);
    fprintf('Carpeta: %s\n', carpeta_base);
    fprintf('============================================\n');

    archivos_hea = dir(fullfile(carpeta_base, '*.hea'));

    if isempty(archivos_hea)
        fprintf('No se encontraron archivos .hea en %s\n', carpeta_base);
        continue
    end

    nombres_registros = cell(numel(archivos_hea),1);

    for r = 1:numel(archivos_hea)
        [~, nombre_sin_ext, ~] = fileparts(archivos_hea(r).name);
        nombres_registros{r} = nombre_sin_ext;
    end

    nombres_registros = unique(nombres_registros);

    %% 4) RECORRER REGISTROS

    for rr = 1:numel(nombres_registros)

        nombre_registro = nombres_registros{rr};

        fprintf('\nProcesando registro: %s (%s)\n', nombre_registro, nombre_base);

        try

            %% 4.1) COMPROBAR ARCHIVOS

            ruta_hea = fullfile(carpeta_base, [nombre_registro '.hea']);
            ruta_dat = fullfile(carpeta_base, [nombre_registro '.dat']);
            ruta_atr = fullfile(carpeta_base, [nombre_registro '.atr']);

            if exist(ruta_hea,'file') ~= 2
                error('No existe el archivo .hea');
            end

            if exist(ruta_dat,'file') ~= 2
                error('No existe el archivo .dat');
            end

            if exist(ruta_atr,'file') ~= 2
                error('No existe el archivo .atr');
            end

            %% 4.2) LEER SEÑAL

            ruta_registro = fullfile(carpeta_base, nombre_registro);

            [sig, Fs_original] = leer_senal_wfdb_python(ruta_registro);

            if isempty(sig)
                error('Registro vacío');
            end

            x = sig(:,1);
            x = double(x(:));

            %% 4.3) REMUESTREO A 500 Hz

            Fs_original_int = round(Fs_original);
            x = resample(x, Fs_new, Fs_original_int);

            Fs = Fs_new;
            dur_total_resamp = length(x)/Fs;

            %% 4.4) FILTRADO PASO BANDA

            [b_filt,a_filt] = butter(2,[0.5 40]/(Fs/2),'bandpass');
            x_filt = filtfilt(b_filt,a_filt,x);

            %% 4.5) DETECCIÓN DE ARTEFACTOS EN VENTANAS DE 2 s

            N_art = dur_artefacto * Fs;
            num_vent_art = floor(length(x_filt)/N_art);

            if num_vent_art < 1
                error('Registro demasiado corto tras remuestreo');
            end

            ptp_vals = zeros(num_vent_art,1);
            rms_vals = zeros(num_vent_art,1);
            der_vals = zeros(num_vent_art,1);
            frac_extreme_vals = zeros(num_vent_art,1);

            for k = 1:num_vent_art

                ini = (k-1)*N_art + 1;
                fin = k*N_art;

                v = x_filt(ini:fin);

                ptp_vals(k) = max(v) - min(v);
                rms_vals(k) = sqrt(mean(v.^2));
                der_vals(k) = median(abs(diff(v)));
            end

            labels = zeros(num_vent_art,1);
            % Etiquetas de calidad:
            %   0 = ventana válida
            %   1 = desconexión
            %   2 = transición tras desconexión
            %   3 = artefacto

            %% 4.6) DETECTAR DESCONEXIÓN

            ref_ptp_all = median(ptp_vals);
            ref_rms_all = median(rms_vals);

            low_ptp = 0.20 * ref_ptp_all;
            low_rms = 0.20 * ref_rms_all;

            idx_disc = find(ptp_vals < low_ptp | rms_vals < low_rms);
            labels(idx_disc) = 1;

            %% 4.7) REFERENCIAS DE SEÑAL NORMAL

            idx_ok_ref = find(labels == 0);

            if isempty(idx_ok_ref)
                error('Todas las ventanas de artefacto han salido malas');
            end

            ref_rms = median(rms_vals(idx_ok_ref));
            ref_ptp = median(ptp_vals(idx_ok_ref));
            ref_der = median(der_vals(idx_ok_ref));

            low_rms_ok  = 0.5 * ref_rms;
            high_rms_ok = 2.0 * ref_rms;

            low_ptp_ok  = 0.5 * ref_ptp;
            high_ptp_ok = 2.0 * ref_ptp;

            high_der_ok = 2.5 * ref_der;

            %% 4.8) DETECTAR TRANSICIÓN DESPUÉS DE DESCONEXIÓN

            k = 1;

            while k <= num_vent_art

                if labels(k) == 0
                    k = k + 1;
                    continue
                end

                while k <= num_vent_art && labels(k) ~= 0
                    k = k + 1;
                end

                stable_count = 0;
                inicio_racha = -1;
                j = k;

                while j <= num_vent_art

                    is_stable = ...
                        rms_vals(j) >= low_rms_ok  && rms_vals(j) <= high_rms_ok && ...
                        ptp_vals(j) >= low_ptp_ok  && ptp_vals(j) <= high_ptp_ok && ...
                        der_vals(j) <= high_der_ok;

                    if is_stable

                        if stable_count == 0
                            inicio_racha = j;
                        end

                        stable_count = stable_count + 1;

                    else

                        stable_count = 0;
                        inicio_racha = -1;
                    end

                    if stable_count >= 3

                        for jj = inicio_racha:j
                            if labels(jj) == 2
                                labels(jj) = 0;
                            end
                        end

                        break
                    end

                    if labels(j) == 0
                        labels(j) = 2;
                    end

                    j = j + 1;
                end
            end

            %% 4.9) DETECTAR ARTEFACTO TIPO PULSOS / AMPLITUD EXTREMA

            idx_no_disc_trans = find(labels == 0);

            if isempty(idx_no_disc_trans)
                error('No quedan ventanas buenas');
            end

            ref_rms_ok2 = median(rms_vals(idx_no_disc_trans));
            amp_thr = 1.3 * ref_rms_ok2;

            for k = 1:num_vent_art

                ini = (k-1)*N_art + 1;
                fin = k*N_art;

                v = x_filt(ini:fin);

                frac_extreme_vals(k) = sum(abs(v) > amp_thr) / length(v);
            end

            ref_frac = median(frac_extreme_vals(idx_no_disc_trans));
            high_frac = max(0.50, 2.0 * ref_frac);

            idx_art = find(labels == 0 & frac_extreme_vals > high_frac);
            labels(idx_art) = 3;

            %% 4.9.1) RELLENAR HUECOS BUENOS CORTOS ENTRE BLOQUES MALOS

            bad = labels ~= 0;
            max_hueco_bueno = 4; % 4 ventanas = 8 s

            k = 1;

            while k <= num_vent_art

                if bad(k)
                    k = k + 1;
                    continue
                end

                ini_hueco = k;

                while k <= num_vent_art && ~bad(k)
                    k = k + 1;
                end

                fin_hueco = k - 1;
                largo_hueco = fin_hueco - ini_hueco + 1;

                if ini_hueco > 1 && fin_hueco < num_vent_art
                    if bad(ini_hueco - 1) && bad(fin_hueco + 1) && largo_hueco <= max_hueco_bueno
                        labels(ini_hueco:fin_hueco) = 3;
                    end
                end
            end

            %% 4.10) CONVERTIR VENTANAS MALAS DE 2 s EN INTERVALOS MALOS

            intervalos_malos = [];

            for k = 1:num_vent_art

                if labels(k) ~= 0
                    intervalos_malos = [intervalos_malos; ...
                        (k-1)*dur_artefacto, k*dur_artefacto]; %#ok<AGROW>
                end
            end

            intervalos_malos = unir_intervalos(intervalos_malos);

            %% 4.11) CONSTRUIR SEGMENTOS LIMPIOS SEGÚN EL TIPO DE BASE

            guardadas_registro = 0;
            descartadas_R_registro = 0;

            %% ========================================================
            % CASO 1: REGISTROS DE RITMO SINUSAL NORMAL
            %% ========================================================

            if strcmp(nombre_base,'BASE_3_SANOS')

                tipo_registro = 'SANO';
                ritmo_ventana = 'SR';
                carpeta_destino = 'SANO';

                fprintf('Tipo de registro: %s\n', tipo_registro);

                [ann_beats, anntype] = leer_anotaciones_latido_wfdb_python(ruta_registro, 'atr');

                if isempty(ann_beats) || isempty(anntype)
                    error('No se pudieron leer anotaciones beat-to-beat en la base sana');
                end

                ann_beats = double(ann_beats(:));
                anntype = string(anntype(:));

                idx_N = anntype == "N";

                if any(idx_N)
                    ann_beats = ann_beats(idx_N);
                end

                if isempty(ann_beats)
                    error('No quedan latidos tras filtrar anotaciones en sano');
                end

                locs_R_global = round(ann_beats * Fs_new / double(Fs_original));
                locs_R_global = unique(locs_R_global);
                locs_R_global = locs_R_global(locs_R_global >= 1 & locs_R_global <= length(x_filt));

                segmentos_limpios = [0 dur_total_resamp];

                for m = 1:size(intervalos_malos,1)

                    segmentos_limpios = restar_intervalo(segmentos_limpios, intervalos_malos(m,:));

                    if isempty(segmentos_limpios)
                        break
                    end
                end

                if isempty(segmentos_limpios)
                    fprintf('No quedan segmentos limpios en %s\n', nombre_registro);
                    continue
                end

                [ventanas_ok, meta_ok, locsR_ok, descartes_local] = extraer_ventanas_desde_segmentos( ...
                    x_filt, Fs, segmentos_limpios, dur_vent_final, locs_R_global);

                for dd = 1:size(descartes_local,1)

                    descartadas_R_registro = descartadas_R_registro + 1;

                    descartadas_R_global(end+1,:) = { ...
                        nombre_base, nombre_registro, ...
                        descartes_local{dd,1}, descartes_local{dd,2}, ...
                        carpeta_destino, descartes_local{dd,3}}; %#ok<AGROW>
                end

                for vv = 1:size(ventanas_ok,1)

                    ventana = ventanas_ok{vv,1};
                    locs_R = locsR_ok{vv,1};

                    t_ini = meta_ok{vv,1};
                    t_fin = meta_ok{vv,2};
                    s = meta_ok{vv,3};
                    w = meta_ok{vv,4};

                    ID_global = ID_global + 1;
                    guardadas_registro = guardadas_registro + 1;

                    nombre_archivo = sprintf('%s_%s_ID_%06d.mat', ...
                        nombre_base, nombre_registro, ID_global);

                    save(fullfile(carpeta_salida, carpeta_destino, nombre_archivo), ...
                        'ventana','locs_R','Fs','nombre_registro','nombre_base', ...
                        's','w','ID_global','t_ini','t_fin', ...
                        'tipo_registro','ritmo_ventana');

                    resumen_global(end+1,:) = { ...
                        nombre_base, nombre_registro, ID_global, ...
                        t_ini, t_fin, tipo_registro, ritmo_ventana, carpeta_destino}; %#ok<AGROW>
                end

            %% ========================================================
            % CASO 2: REGISTROS CON FA O RITMO SINUSAL ANOTADO
            %% ========================================================

            else

                [ann_ritmo, comments] = leer_anotaciones_wfdb_python(ruta_registro, 'atr');

                if isempty(ann_ritmo) || isempty(comments)
                    error('No hay anotaciones de ritmo');
                end

                ann_qrs = leer_anotaciones_muestra_wfdb_python(ruta_registro, 'qrs');

                if isempty(ann_qrs)
                    error('No se pudieron leer anotaciones .qrs');
                end

                locs_R_global = round(double(ann_qrs(:)) * Fs_new / double(Fs_original));
                locs_R_global = unique(locs_R_global);
                locs_R_global = locs_R_global(locs_R_global >= 1 & locs_R_global <= length(x_filt));

                ritmos_validos = {'(AFIB','(N','(NSR','(SR'};
                tiempos_ritmo = [];
                ritmos = {};

                for i = 1:min(length(ann_ritmo), length(comments))

                    txt = strtrim(comments{i});

                    if any(strcmp(txt, ritmos_validos))
                        tiempos_ritmo(end+1,1) = ann_ritmo(i) / double(Fs_original); %#ok<AGROW>
                        ritmos{end+1,1} = txt; %#ok<AGROW>
                    end
                end

                if isempty(tiempos_ritmo)
                    error('No se encontraron etiquetas de ritmo útiles');
                end

                [tiempos_ritmo, idx_sort] = sort(tiempos_ritmo);
                ritmos = ritmos(idx_sort);

                hayAF = any(strcmp(ritmos,'(AFIB'));
                haySR = any(strcmp(ritmos,'(N')) || ...
                        any(strcmp(ritmos,'(NSR')) || ...
                        any(strcmp(ritmos,'(SR'));

                if hayAF && haySR
                    tipo_registro = 'FA_PAROXISTICA';
                elseif hayAF
                    tipo_registro = 'FA_PERSISTENTE';
                elseif haySR
                    tipo_registro = 'SANO';
                else
                    tipo_registro = 'REVISAR';
                end

                fprintf('Tipo de registro: %s\n', tipo_registro);

                for i = 1:length(tiempos_ritmo)

                    seg_ini = tiempos_ritmo(i);

                    if i < length(tiempos_ritmo)
                        seg_fin = tiempos_ritmo(i+1);
                    else
                        seg_fin = dur_total_resamp;
                    end

                    if seg_fin <= seg_ini
                        continue
                    end

                    ritmo_ventana = mapear_ritmo(ritmos{i});

                    if strcmp(ritmo_ventana,'OTHER')
                        continue
                    end

                    carpeta_destino = decidir_carpeta(tipo_registro, ritmo_ventana);

                    if isempty(carpeta_destino)
                        continue
                    end

                    segmentos_limpios = [seg_ini seg_fin];

                    for m = 1:size(intervalos_malos,1)

                        segmentos_limpios = restar_intervalo(segmentos_limpios, intervalos_malos(m,:));

                        if isempty(segmentos_limpios)
                            break
                        end
                    end

                    if isempty(segmentos_limpios)
                        continue
                    end

                    [ventanas_ok, meta_ok, locsR_ok, descartes_local] = extraer_ventanas_desde_segmentos( ...
                        x_filt, Fs, segmentos_limpios, dur_vent_final, locs_R_global);

                    for dd = 1:size(descartes_local,1)

                        descartadas_R_registro = descartadas_R_registro + 1;

                        descartadas_R_global(end+1,:) = { ...
                            nombre_base, nombre_registro, ...
                            descartes_local{dd,1}, descartes_local{dd,2}, ...
                            carpeta_destino, descartes_local{dd,3}}; %#ok<AGROW>
                    end

                    for vv = 1:size(ventanas_ok,1)

                        ventana = ventanas_ok{vv,1};
                        locs_R = locsR_ok{vv,1};

                        t_ini = meta_ok{vv,1};
                        t_fin = meta_ok{vv,2};
                        s = meta_ok{vv,3};
                        w = meta_ok{vv,4};

                        ID_global = ID_global + 1;
                        guardadas_registro = guardadas_registro + 1;

                        nombre_archivo = sprintf('%s_%s_ID_%06d.mat', ...
                            nombre_base, nombre_registro, ID_global);

                        save(fullfile(carpeta_salida, carpeta_destino, nombre_archivo), ...
                            'ventana','locs_R','Fs','nombre_registro','nombre_base', ...
                            's','w','ID_global','t_ini','t_fin', ...
                            'tipo_registro','ritmo_ventana');

                        resumen_global(end+1,:) = { ...
                            nombre_base, nombre_registro, ID_global, ...
                            t_ini, t_fin, tipo_registro, ritmo_ventana, carpeta_destino}; %#ok<AGROW>
                    end
                end
            end

            %% 4.12) RESUMEN POR REGISTRO

            fprintf('Ventanas artefacto buenas       : %d\n', sum(labels == 0));
            fprintf('Desconexion                     : %d\n', sum(labels == 1));
            fprintf('Transicion                      : %d\n', sum(labels == 2));
            fprintf('Artefacto                       : %d\n', sum(labels == 3));
            fprintf('Descartadas por calidad/ritmo   : %d\n', descartadas_R_registro);
            fprintf('Ventanas finales guardadas      : %d\n', guardadas_registro);

        catch ME

            fprintf('ERROR en %s (%s): %s\n', nombre_registro, nombre_base, ME.message);

            if ~isempty(ME.stack)
                fprintf('Linea aproximada: %d\n', ME.stack(1).line);
            end

            errores_global(end+1,:) = {nombre_base, nombre_registro, ME.message}; %#ok<AGROW>

            continue
        end
    end
end

%% 5) GUARDAR RESUMEN Y ERRORES

if ~isempty(resumen_global)

    T_resumen = cell2table(resumen_global, 'VariableNames', ...
        {'Base','Registro','ID_Ventana','TiempoInicio','TiempoFin', ...
         'TipoRegistro','RitmoVentana','CarpetaDestino'});

    writetable(T_resumen, fullfile(carpeta_salida, 'resumen_ventanas_global.csv'));
end

if ~isempty(errores_global)

    T_errores = cell2table(errores_global, 'VariableNames', ...
        {'Base','Registro','MensajeError'});

    writetable(T_errores, fullfile(carpeta_salida, 'errores_procesamiento_global.csv'));
end

if ~isempty(descartadas_R_global)

    T_desc_R = cell2table(descartadas_R_global, 'VariableNames', ...
        {'Base','Registro','TiempoInicio','TiempoFin','CarpetaDestino','Motivo'});

    writetable(T_desc_R, fullfile(carpeta_salida, 'ventanas_descartadas_global.csv'));
end

fprintf('\nProceso global terminado.\n');
fprintf('Último ID global: %d\n', ID_global);

%% ============================================================
% FUNCIÓN LOCAL: EXTRAER VENTANAS DESDE SEGMENTOS LIMPIOS
%% ============================================================

function [ventanas_ok, meta_ok, locsR_ok, descartes] = extraer_ventanas_desde_segmentos( ...
    sig, Fs, segmentos_limpios, dur_vent, locs_R_global)

ventanas_ok = {};
meta_ok = {};
locsR_ok = {};
descartes = {};

N_final = dur_vent * Fs;

for s = 1:size(segmentos_limpios,1)

    limpio_ini = segmentos_limpios(s,1);
    limpio_fin = segmentos_limpios(s,2);

    dur_limpio = limpio_fin - limpio_ini;

    if dur_limpio < dur_vent
        continue
    end

    num_subvent = floor(dur_limpio / dur_vent);

    for w = 1:num_subvent

        t_ini = limpio_ini + (w-1)*dur_vent;
        t_fin = t_ini + dur_vent;

        ini = round(t_ini * Fs) + 1;
        fin = ini + N_final - 1;

        if ini < 1 || fin > length(sig)
            continue
        end

        ventana = sig(ini:fin);

        %% R anotados dentro de la ventana en coordenadas locales

        if nargin >= 5 && ~isempty(locs_R_global)

            idx_R = locs_R_global >= ini & locs_R_global <= fin;

            locs_R_anot_abs = locs_R_global(idx_R);
            locs_R_anot = locs_R_anot_abs - ini + 1;

        else
            locs_R_anot = [];
        end

        %% Método híbrido conservador

        [locs_R, motivo_R] = DetectarPICOSR(ventana, locs_R_anot, Fs);

        if isempty(locs_R) || numel(locs_R) < 3

            if isempty(motivo_R)
                motivo_R = 'Menos de 3 picos R tras método híbrido';
            end

            descartes(end+1,:) = {t_ini, t_fin, motivo_R}; %#ok<AGROW>
            continue
        end

        %% Control RR

        RR = diff(locs_R) / Fs;

        if isempty(RR) || numel(RR) < 2
            descartes(end+1,:) = {t_ini, t_fin, 'RR insuficientes'}; %#ok<AGROW>
            continue
        end

        RR_validos = RR(isfinite(RR));

        if isempty(RR_validos)
            descartes(end+1,:) = {t_ini, t_fin, 'RR no finitos'}; %#ok<AGROW>
            continue
        end

        HR_mean = 60 / mean(RR_validos, 'omitnan');

        if ~isfinite(HR_mean) || HR_mean < 25 || HR_mean > 250
            descartes(end+1,:) = {t_ini, t_fin, 'FC fuera de rango'}; %#ok<AGROW>
            continue
        end

        ventanas_ok(end+1,1) = {ventana}; %#ok<AGROW>
        meta_ok(end+1,:) = {t_ini, t_fin, s, w}; %#ok<AGROW>
        locsR_ok(end+1,1) = {locs_R}; %#ok<AGROW>
    end
end

end

function ultimo_id = buscar_ultimo_id_global(carpeta_salida)
    subcarpetas = {'SANO','FA_PERSISTENTE','FA_PAROXISTICA_RS','FA_PAROXISTICA_FA'};
    ultimo_id = 0;

    for ii = 1:numel(subcarpetas)
        carpeta_actual = fullfile(carpeta_salida, subcarpetas{ii});
        archivos = dir(fullfile(carpeta_actual, '*.mat'));

        for k = 1:numel(archivos)
            nombre = archivos(k).name;
            token = regexp(nombre, 'ID_(\d+)\.mat$', 'tokens', 'once');
            if ~isempty(token)
                id_actual = str2double(token{1});
                if id_actual > ultimo_id
                    ultimo_id = id_actual;
                end
            end
        end
    end
end

function [sig, Fs] = leer_senal_wfdb_python(ruta_registro)

    wfdb = py.importlib.import_module('wfdb');

    % Leer solo el canal 0
    out = wfdb.rdsamp(ruta_registro, pyargs('channels', py.list({int32(0)})));

    sig_py = out{1};      % señal
    campos = out{2};      % metadatos

    Fs = [];
    try
        if isa(campos, 'py.dict')
            Fs = double(campos{'fs'});
        else
            Fs = double(campos.fs);
        end
    catch ME
        error('No se pudo obtener la frecuencia de muestreo: %s', ME.message);
    end
    sig = py_numpy_array_to_mat(sig_py);

    if isempty(sig)
        error('No se pudo convertir p_signal a matriz MATLAB.');
    end
end

function [ann, anntype] = leer_anotaciones_latido_wfdb_python(ruta_registro, anotador)

    wfdb = py.importlib.import_module('wfdb');
    ann_py = wfdb.rdann(ruta_registro, anotador);

    np = py.importlib.import_module('numpy');
    sample_np = np.asarray(ann_py.sample, pyargs('dtype', np.int64));
    ann = py_numpy_array_to_mat(sample_np);
    ann = ann(:);

    symbol_list = cell(ann_py.symbol);

    n = min(numel(ann), numel(symbol_list));  
    ann = ann(1:n);

    anntype = cell(n,1);
    for i = 1:n
        elem = symbol_list{i};
        if isempty(elem) || strcmp(class(elem),'py.NoneType')
            anntype{i} = '';
        else
            anntype{i} = char(string(elem));
        end
    end

end

function [ann, comments] = leer_anotaciones_wfdb_python(ruta_registro, extension_anot)
    np = py.importlib.import_module('numpy');
    an = py.wfdb.rdann(ruta_registro, extension_anot);

    sample_np = np.asarray(an.sample, pyargs('dtype', np.int64));
    ann = py_numpy_array_to_mat(sample_np);
    ann = ann(:);

    aux_list = cell(an.aux_note);
    comments = cell(numel(aux_list),1);

    for ii = 1:numel(aux_list)
        elem = aux_list{ii};
        if isempty(elem) || strcmp(class(elem), 'py.NoneType')
            comments{ii} = '';
        else
            comments{ii} = strtrim(char(string(elem)));
        end
    end
end


function ann = leer_anotaciones_muestra_wfdb_python(ruta_registro, anotador)

ann = [];

try
    wfdb = py.importlib.import_module('wfdb');
    ann_py = wfdb.rdann(ruta_registro, anotador);

    np = py.importlib.import_module('numpy');

    sample_np = np.asarray(ann_py.sample, pyargs('dtype', np.int64));

    ann = py_numpy_array_to_mat(sample_np);
    ann = ann(:);

catch ME
    warning('Error leyendo muestras WFDB (%s): %s', anotador, ME.message);
    ann = [];
end

end

function A = py_numpy_array_to_mat(py_array)

    np = py.importlib.import_module('numpy');

    arr = np.asarray(py_array, pyargs('dtype', np.float64, 'order', 'C'));
    shape = cell(arr.shape);
    dims = cellfun(@double, shape);

    try
        v = double(py.array.array('d', arr.flatten().tolist()));
    catch
        data_cell = cell(arr.flatten().tolist());
        v = cellfun(@double, data_cell);
    end

    if isempty(dims)
        A = [];
    elseif numel(dims) == 1
        A = reshape(v, [dims(1), 1]);
    elseif numel(dims) == 2
        A = reshape(v, fliplr(dims))';
    else
        error('Dimensionalidad no soportada en py_numpy_array_to_mat.');
    end
end


function intervalos_unidos = unir_intervalos(intervalos)
    if isempty(intervalos)
        intervalos_unidos = [];
        return
    end

    intervalos = sortrows(intervalos,1);
    intervalos_unidos = intervalos(1,:);

    for m = 2:size(intervalos,1)
        actual_ini = intervalos(m,1);
        actual_fin = intervalos(m,2);
        ultimo_fin = intervalos_unidos(end,2);

        if actual_ini <= ultimo_fin
            intervalos_unidos(end,2) = max(ultimo_fin, actual_fin);
        else
            intervalos_unidos = [intervalos_unidos; actual_ini actual_fin];
        end
    end
end


function segmentos_out = restar_intervalo(segmentos_in, intervalo_malo)
    malo_ini = intervalo_malo(1);
    malo_fin = intervalo_malo(2);

    segmentos_out = [];

    for s = 1:size(segmentos_in,1)
        limpio_ini = segmentos_in(s,1);
        limpio_fin = segmentos_in(s,2);

        if malo_fin <= limpio_ini || malo_ini >= limpio_fin
            segmentos_out = [segmentos_out; limpio_ini limpio_fin];

        elseif malo_ini <= limpio_ini && malo_fin >= limpio_fin
            % desaparece entero

        elseif malo_ini <= limpio_ini && malo_fin < limpio_fin
            segmentos_out = [segmentos_out; malo_fin limpio_fin];

        elseif malo_ini > limpio_ini && malo_fin >= limpio_fin
            segmentos_out = [segmentos_out; limpio_ini malo_ini];

        elseif malo_ini > limpio_ini && malo_fin < limpio_fin
            segmentos_out = [segmentos_out; limpio_ini malo_ini];
            segmentos_out = [segmentos_out; malo_fin limpio_fin];
        end
    end
end

function ritmo = mapear_ritmo(txt)
    if strcmp(txt,'(AFIB')
        ritmo = 'FA';
    elseif strcmp(txt,'(N') || strcmp(txt,'(NSR') || strcmp(txt,'(SR')
        ritmo = 'RS';
    else
        ritmo = 'OTHER';
    end
end


function carpeta = decidir_carpeta(tipo_registro, ritmo_intervalo)
    carpeta = '';

    if strcmp(tipo_registro,'SANO')
        if strcmp(ritmo_intervalo,'RS')
            carpeta = 'SANO';
        end
    elseif strcmp(tipo_registro,'FA_PERSISTENTE')
        if strcmp(ritmo_intervalo,'FA')
            carpeta = 'FA_PERSISTENTE';
        end
    elseif strcmp(tipo_registro,'FA_PAROXISTICA')
        if strcmp(ritmo_intervalo,'RS')
            carpeta = 'FA_PAROXISTICA_RS';
        elseif strcmp(ritmo_intervalo,'FA')
            carpeta = 'FA_PAROXISTICA_FA';
        end
    end
end
