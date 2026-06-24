clear
clc
close all

%% ============================================================
% ANÁLISIS DE TRANSICIONES RS -> FA
%
% Este script analiza transiciones de ritmo sinusal a fibrilación auricular
% en registros de larga duración. Para cada transición válida se estudia la
% frecuencia dominante (DF) del residual auricular en cuatro ventanas de 5 s
% previas al inicio del episodio de FA.
%
% Estrategia general:
%   - Se leen la señal ECG y las anotaciones de ritmo de los registros.
%   - Se descarta el registro 112.
%   - Se detectan y excluyen ventanas con artefactos.
%   - Se utilizan las anotaciones .qrs como referencia principal para los
%     picos R.
%   - DetectarPICOSR complementa los picos R cuando existen huecos largos.
%   - Se realiza un ajuste final de los picos R al máximo local.
%   - Se cancela el complejo QRS-T mediante plantilla mediana estrecha,
%     con el objetivo de preservar la onda P.
%   - Se calcula la DF del residual auricular en ventanas previas al inicio
%     de la FA.
%
% Ventanas analizadas respecto al inicio de la FA:
%   - Ventana alejada 180 s: -180 a -175 s
%   - Ventana alejada 60 s:   -60 a  -55 s
%   - Ventana penúltima:      -10 a   -5 s
%   - Ventana final:           -5 a    0 s
%
% Comparaciones realizadas:
%   - Delta cercano = DF(-5 a 0 s) - DF(-10 a -5 s)
%   - Delta 60 s    = DF(-5 a 0 s) - DF(-60 a -55 s)
%   - Delta 180 s   = DF(-5 a 0 s) - DF(-180 a -175 s)
%
% Salidas:
%   - resultados_RS_FA_ventanas.xlsx
%   - resumen_estadistico_RS_FA_ventanas.xlsx
%   - boxplots de las ventanas analizadas
%   - gráficos pareados
%   - histogramas de los deltas de DF
%
% Funciones locales incluidas:
%   - extraer_ventana_relativa_fin
%   - leer_senal_wfdb_python
%   - leer_anotaciones_wfdb_python
%   - leer_anotaciones_muestra_wfdb_python
%   - py_numpy_array_to_mat
%   - unir_intervalos
%   - mapear_ritmo
%   - limpiar_locs_local
%   - cancelar_QRST_plantilla_medianaRS
%   - frecuencia_dominante2_RS
%
% Función externa necesaria:
%   - DetectarPICOSR.m
%
% Requisitos:
%   - Python con wfdb y numpy.
%   - Signal Processing Toolbox.
%   - Statistics and Machine Learning Toolbox.
%% ============================================================

%% CONFIGURACIÓN GENERAL

bases = { ...
    struct('carpeta','C:\Users\Emma\Documents\MATLAB\long-term-af-database-1.0.0\files', ...
           'nombre','BASE_1'), ...
    struct('carpeta','C:\Users\Emma\Documents\MATLAB\mit-bih-atrial-fibrillation-database-1.0.0\files', ...
           'nombre','BASE_2') ...
};

carpeta_out = 'C:\Users\Emma\Documents\MATLAB\Petrutiu_RS_FA_sinartefactosDEFINITIVO';

if ~exist(carpeta_out, 'dir')
    mkdir(carpeta_out);
end

Fs_new = 500;

dur_cancel_RS = 180;   % segundos antes del inicio de FA para cancelar QRS-T

% Ventanas que se van a analizar respecto al inicio de FA
vent_alejada180_ini = -180;
vent_alejada180_fin = -175;

vent_alejada60_ini = -60;
vent_alejada60_fin = -55;

vent_pen_ini = -10;
vent_pen_fin = -5;

vent_ult_ini = -5;
vent_ult_fin = 0;

min_R_RS_cancel = 15;

registros_excluir = {'112'};

resultados_RS_FA = {};

%% ============================================================
% RECORRER BASES
%% ============================================================

for bb = 1:numel(bases)

    carpeta_base = bases{bb}.carpeta;
    nombre_base  = bases{bb}.nombre;

    archivos_hea = dir(fullfile(carpeta_base, '*.hea'));

    fprintf('\n============================================\n');
    fprintf('Procesando %s\n', nombre_base);
    fprintf('Carpeta: %s\n', carpeta_base);
    fprintf('Registros encontrados: %d\n', numel(archivos_hea));
    fprintf('============================================\n');

    for rr = 1:numel(archivos_hea)

        [~, nombre_registro, ~] = fileparts(archivos_hea(rr).name);

        if any(strcmp(nombre_registro, registros_excluir))
            fprintf('\nRegistro %s excluido del análisis.\n', nombre_registro);
            continue
        end

        ruta_registro = fullfile(carpeta_base, nombre_registro);

        fprintf('\nProcesando %s | %s\n', nombre_base, nombre_registro);

        try

            %% ============================================================
            % 1) LEER SEÑAL
            %% ============================================================

            [sig, Fs_original] = leer_senal_wfdb_python(ruta_registro);

            if isempty(sig)
                fprintf('  Señal vacía. Se omite.\n');
                continue
            end

            x = sig(:,1);
            x = double(x(:));

            %% ============================================================
            % 2) REMUESTREO A 500 Hz
            %% ============================================================

            x = resample(x, Fs_new, round(Fs_original));
            Fs = Fs_new;

            %% ============================================================
            % 2.5) DETECCIÓN DE ARTEFACTOS
            %
            % Se utiliza un filtrado independiente entre 0,5 y 40 Hz
            % exclusivamente para evaluar la calidad de la señal.
            % El filtrado 0,5-20 Hz del análisis de RS se mantiene aparte.
            %% ============================================================

            [b_calidad, a_calidad] = butter(2, [0.5 40] / (Fs/2), 'bandpass');
            x_calidad = filtfilt(b_calidad, a_calidad, x);

            dur_artefacto = 2;
            N_art = dur_artefacto * Fs;
            num_vent_art = floor(length(x_calidad) / N_art);

            if num_vent_art < 1
                fprintf('  Registro demasiado corto para evaluar artefactos. Se omite.\n');
                continue
            end

            ptp_vals = zeros(num_vent_art,1);
            rms_vals = zeros(num_vent_art,1);
            der_vals = zeros(num_vent_art,1);
            frac_extreme_vals = zeros(num_vent_art,1);

            for k_art = 1:num_vent_art

                ini_art = (k_art-1)*N_art + 1;
                fin_art = k_art*N_art;

                v_art = x_calidad(ini_art:fin_art);

                ptp_vals(k_art) = max(v_art) - min(v_art);
                rms_vals(k_art) = sqrt(mean(v_art.^2));
                der_vals(k_art) = median(abs(diff(v_art)));
            end

            labels = zeros(num_vent_art,1);
            % 0 buena, 1 desconexión, 2 transición, 3 artefacto

            %% Detectar desconexión

            ref_ptp_all = median(ptp_vals);
            ref_rms_all = median(rms_vals);

            low_ptp = 0.20 * ref_ptp_all;
            low_rms = 0.20 * ref_rms_all;

            idx_disc = find(ptp_vals < low_ptp | rms_vals < low_rms);
            labels(idx_disc) = 1;

            %% Referencias de señal normal

            idx_ok_ref = find(labels == 0);

            if isempty(idx_ok_ref)
                fprintf('  Todas las ventanas de calidad han salido malas. Se omite.\n');
                continue
            end

            ref_rms = median(rms_vals(idx_ok_ref));
            ref_ptp = median(ptp_vals(idx_ok_ref));
            ref_der = median(der_vals(idx_ok_ref));

            low_rms_ok  = 0.5 * ref_rms;
            high_rms_ok = 2.0 * ref_rms;

            low_ptp_ok  = 0.5 * ref_ptp;
            high_ptp_ok = 2.0 * ref_ptp;

            high_der_ok = 2.5 * ref_der;

            %% Detectar transición después de desconexión

            k_art = 1;

            while k_art <= num_vent_art

                if labels(k_art) == 0
                    k_art = k_art + 1;
                    continue
                end

                while k_art <= num_vent_art && labels(k_art) ~= 0
                    k_art = k_art + 1;
                end

                stable_count = 0;
                inicio_racha = -1;
                j_art = k_art;

                while j_art <= num_vent_art

                    is_stable = ...
                        rms_vals(j_art) >= low_rms_ok && ...
                        rms_vals(j_art) <= high_rms_ok && ...
                        ptp_vals(j_art) >= low_ptp_ok && ...
                        ptp_vals(j_art) <= high_ptp_ok && ...
                        der_vals(j_art) <= high_der_ok;

                    if is_stable

                        if stable_count == 0
                            inicio_racha = j_art;
                        end

                        stable_count = stable_count + 1;

                    else

                        stable_count = 0;
                        inicio_racha = -1;
                    end

                    if stable_count >= 3

                        for jj_art = inicio_racha:j_art
                            if labels(jj_art) == 2
                                labels(jj_art) = 0;
                            end
                        end

                        break
                    end

                    if labels(j_art) == 0
                        labels(j_art) = 2;
                    end

                    j_art = j_art + 1;
                end
            end

            %% Detectar pulsos o amplitud extrema

            idx_no_disc_trans = find(labels == 0);

            if isempty(idx_no_disc_trans)
                fprintf('  No quedan ventanas buenas tras evaluar desconexiones. Se omite.\n');
                continue
            end

            ref_rms_ok2 = median(rms_vals(idx_no_disc_trans));
            amp_thr = 1.3 * ref_rms_ok2;

            for k_art = 1:num_vent_art

                ini_art = (k_art-1)*N_art + 1;
                fin_art = k_art*N_art;

                v_art = x_calidad(ini_art:fin_art);

                frac_extreme_vals(k_art) = ...
                    sum(abs(v_art) > amp_thr) / length(v_art);
            end

            ref_frac = median(frac_extreme_vals(idx_no_disc_trans));
            high_frac = max(0.50, 2.0 * ref_frac);

            idx_art = find(labels == 0 & frac_extreme_vals > high_frac);
            labels(idx_art) = 3;

            %% Rellenar huecos buenos cortos entre bloques malos

            bad = labels ~= 0;
            max_hueco_bueno = 4; % 4 ventanas = 8 s

            k_art = 1;

            while k_art <= num_vent_art

                if bad(k_art)
                    k_art = k_art + 1;
                    continue
                end

                ini_hueco = k_art;

                while k_art <= num_vent_art && ~bad(k_art)
                    k_art = k_art + 1;
                end

                fin_hueco = k_art - 1;
                largo_hueco = fin_hueco - ini_hueco + 1;

                if ini_hueco > 1 && fin_hueco < num_vent_art

                    if bad(ini_hueco - 1) && ...
                            bad(fin_hueco + 1) && ...
                            largo_hueco <= max_hueco_bueno

                        labels(ini_hueco:fin_hueco) = 3;
                    end
                end
            end

            %% Convertir ventanas malas de 2 s en intervalos malos

            intervalos_malos = [];

            for k_art = 1:num_vent_art

                if labels(k_art) ~= 0

                    intervalos_malos = [intervalos_malos; ...
                        (k_art-1)*dur_artefacto, ...
                        k_art*dur_artefacto]; %#ok<AGROW>
                end
            end

            intervalos_malos = unir_intervalos(intervalos_malos);

            %% ============================================================
            % 3) FILTRADO PARA ANÁLISIS DE RS
            %% ============================================================

            [b_RS, a_RS] = butter(2, [0.5 20] / (Fs/2), 'bandpass');
            x_filt_RS = filtfilt(b_RS, a_RS, x);

            %% ============================================================
            % 4) ANOTACIONES DE RITMO
            %% ============================================================

            [ann_ritmo, comments] = leer_anotaciones_wfdb_python(ruta_registro, 'atr');

            if isempty(ann_ritmo) || isempty(comments)
                fprintf('  Sin anotaciones de ritmo .atr. Se omite.\n');
                continue
            end

            %% ============================================================
            % 5) ANOTACIONES .qrs GLOBALES
            %% ============================================================

            ann_qrs = leer_anotaciones_muestra_wfdb_python(ruta_registro, 'qrs');

            if isempty(ann_qrs)
                fprintf('  Sin anotaciones .qrs. Se omite.\n');
                continue
            end

            locs_R_qrs_global = round(double(ann_qrs(:)) * Fs_new / double(Fs_original));
            locs_R_qrs_global = unique(locs_R_qrs_global);
            locs_R_qrs_global = locs_R_qrs_global(locs_R_qrs_global >= 1 & locs_R_qrs_global <= length(x));

            if numel(locs_R_qrs_global) < 3
                fprintf('  Pocos R .qrs globales. Se omite.\n');
                continue
            end

            %% ============================================================
            % 6) PROCESAR ANOTACIONES DE RITMO
            %% ============================================================

            % IMPORTANTE:
            % Estas son etiquetas originales de PhysioNet. No cambiarlas.
            %
            % Las anotaciones consecutivas equivalentes se unifican después
            % de ordenarlas, conservando la primera de cada episodio.
            ritmos_validos = {'(AFIB','(N','(NSR','(SR'};

            tiempos_ritmo = [];
            ritmos = {};

            for i = 1:min(length(ann_ritmo), length(comments))

                txt = strtrim(comments{i});

                if any(strcmp(txt, ritmos_validos))
                    tiempos_ritmo(end+1,1) = ann_ritmo(i) / double(Fs_original); %#ok<SAGROW>
                    ritmos{end+1,1} = txt; %#ok<SAGROW>
                end
            end

            if isempty(tiempos_ritmo)
                fprintf('  No hay ritmos válidos. Se omite.\n');
                continue
            end

            [tiempos_ritmo, idx] = sort(tiempos_ritmo);
            ritmos = ritmos(idx);

            %% ============================================================
            % 6.1) UNIFICAR ANOTACIONES CONSECUTIVAS DEL MISMO RITMO
            %
            % Las anotaciones consecutivas que representan el mismo ritmo
            % se consideran parte de un único episodio continuo.
            %
            % Ejemplo:
            %   0 s    (N
            %   400 s  (NSR
            %   500 s  (AFIB
            %
            % Se transforma en:
            %   0 s    (N
            %   500 s  (AFIB
            %
            % Así se conserva el inicio real del episodio continuo de RS.
            %% ============================================================

            tiempos_ritmo_limpios = [];
            ritmos_limpios = {};
            n_repetidas_eliminadas = 0;

            for i = 1:numel(ritmos)

                ritmo_actual_mapeado = mapear_ritmo(ritmos{i});

                if isempty(ritmos_limpios)

                    tiempos_ritmo_limpios(end+1,1) = ...
                        tiempos_ritmo(i); %#ok<AGROW>

                    ritmos_limpios{end+1,1} = ...
                        ritmos{i}; %#ok<AGROW>

                else

                    ritmo_anterior_mapeado = ...
                        mapear_ritmo(ritmos_limpios{end});

                    if strcmp(ritmo_actual_mapeado, ...
                            ritmo_anterior_mapeado)

                        n_repetidas_eliminadas = ...
                            n_repetidas_eliminadas + 1;

                        fprintf(['  Anotacion repetida unificada: ' ...
                            '%s en %.1f s\n'], ...
                            ritmos{i}, ...
                            tiempos_ritmo(i));

                    else

                        tiempos_ritmo_limpios(end+1,1) = ...
                            tiempos_ritmo(i); %#ok<AGROW>

                        ritmos_limpios{end+1,1} = ...
                            ritmos{i}; %#ok<AGROW>
                    end
                end
            end

            tiempos_ritmo = tiempos_ritmo_limpios;
            ritmos = ritmos_limpios;

            fprintf('  Anotaciones repetidas unificadas: %d\n', ...
                n_repetidas_eliminadas);

            fprintf('  Episodios de ritmo tras la unificacion: %d\n', ...
                numel(ritmos));

            if numel(tiempos_ritmo) < 2

                fprintf(['  No quedan cambios de ritmo suficientes ' ...
                    'tras unificar repeticiones.\n']);

                continue
            end

            %% ============================================================
            % 7) ANÁLISIS RS -> FA
            %% ============================================================

            n_RS_FA_reg = 0;

            for i = 1:length(tiempos_ritmo)-1

                ritmo_actual = mapear_ritmo(ritmos{i});
                ritmo_sig    = mapear_ritmo(ritmos{i+1});

                % Compatibilidad por si mapear_ritmo devuelve AF/SR o FA/RS
                es_RS_actual = strcmp(ritmo_actual,'RS') || strcmp(ritmo_actual,'SR');
                es_FA_sig    = strcmp(ritmo_sig,'FA')    || strcmp(ritmo_sig,'AF');

                if ~(es_RS_actual && es_FA_sig)
                    continue
                end

                t_ini_RS = tiempos_ritmo(i);
                t_ini_FA = tiempos_ritmo(i+1);

                if (t_ini_FA - t_ini_RS) < dur_cancel_RS
                    continue
                end

                %% --------------------------------------------------------
                % Ventana local de 180 s antes del inicio de FA
                %% --------------------------------------------------------

                t0_local = t_ini_FA - dur_cancel_RS;
                t1_local = t_ini_FA;

                ini_local = round(t0_local * Fs) + 1;
                fin_local = ini_local + dur_cancel_RS * Fs - 1;

                if ini_local < 1 || fin_local > length(x_filt_RS)
                    continue
                end

                %% --------------------------------------------------------
                % Comprobar artefactos en las cuatro ventanas analizadas
                %
                % Solo se descarta la transición si un intervalo malo se
                % solapa con:
                %   -180 a -175 s
                %    -60 a  -55 s
                %    -10 a   -5 s
                %     -5 a    0 s
                %% --------------------------------------------------------

                ventanas_calidad = [ ...
                    t_ini_FA + vent_alejada180_ini, ...
                    t_ini_FA + vent_alejada180_fin; ...
                    t_ini_FA + vent_alejada60_ini, ...
                    t_ini_FA + vent_alejada60_fin; ...
                    t_ini_FA + vent_pen_ini, ...
                    t_ini_FA + vent_pen_fin; ...
                    t_ini_FA + vent_ult_ini, ...
                    t_ini_FA + vent_ult_fin];

                hay_artefacto_ventanas = false;

                for vv_art = 1:size(ventanas_calidad,1)

                    t_ini_comprobar = ventanas_calidad(vv_art,1);
                    t_fin_comprobar = ventanas_calidad(vv_art,2);

                    for aa_art = 1:size(intervalos_malos,1)

                        if intervalos_malos(aa_art,1) < t_fin_comprobar && ...
                                intervalos_malos(aa_art,2) > t_ini_comprobar

                            hay_artefacto_ventanas = true;
                            break
                        end
                    end

                    if hay_artefacto_ventanas
                        break
                    end
                end

                if hay_artefacto_ventanas

                    fprintf(['  RS->FA omitida t=%.2f s: alguna ventana ' ...
                        'analizada contiene artefactos.\n'], t_ini_FA);

                    continue
                end

                ecg_local_RS = x_filt_RS(ini_local:fin_local);
                ecg_local_RS = ecg_local_RS(:);

                if length(ecg_local_RS) < dur_cancel_RS * Fs
                    continue
                end

                %% --------------------------------------------------------
                % R .qrs locales
                %% --------------------------------------------------------

                idx_qrs_local = locs_R_qrs_global >= ini_local & locs_R_qrs_global <= fin_local;
                locs_qrs_local = locs_R_qrs_global(idx_qrs_local) - ini_local + 1;
                locs_qrs_local = limpiar_locs_local(locs_qrs_local, length(ecg_local_RS));

                %% --------------------------------------------------------
                % R híbridos locales
                %% --------------------------------------------------------

                [locs_R_local, motivo_R] = DetectarPICOSR(ecg_local_RS, locs_qrs_local, Fs);

                if isempty(locs_R_local) || numel(locs_R_local) < min_R_RS_cancel
                    fprintf('  RS->FA omitida t=%.2f s: R insuficientes. %s\n', ...
                        t_ini_FA, motivo_R);
                    continue
                end

                %% --------------------------------------------------------
                % Cancelar QRS-T localmente preservando P
                %% --------------------------------------------------------

                atrial_local_RS = cancelar_QRST_plantilla_medianaRS(ecg_local_RS, locs_R_local, Fs);

                if isempty(atrial_local_RS) || any(~isfinite(atrial_local_RS))
                    fprintf('  RS->FA omitida t=%.2f s: residual RS no válido.\n', t_ini_FA);
                    continue
                end

                atrial_local_RS = atrial_local_RS(:);

                %% --------------------------------------------------------
                % Extraer las cuatro ventanas respecto al inicio de FA
                %
                % Como atrial_local_RS va desde:
                %   t_ini_FA - 180 s  hasta  t_ini_FA
                %
                % El instante relativo 0 corresponde al final de la ventana.
                %% --------------------------------------------------------

                [seg_alejado180, ok_alejado180] = extraer_ventana_relativa_fin( ...
                    atrial_local_RS, Fs, vent_alejada180_ini, vent_alejada180_fin);

                [seg_alejado60, ok_alejado60] = extraer_ventana_relativa_fin( ...
                    atrial_local_RS, Fs, vent_alejada60_ini, vent_alejada60_fin);

                [seg_penult, ok_pen] = extraer_ventana_relativa_fin( ...
                    atrial_local_RS, Fs, vent_pen_ini, vent_pen_fin);

                [seg_ult, ok_ult] = extraer_ventana_relativa_fin( ...
                    atrial_local_RS, Fs, vent_ult_ini, vent_ult_fin);

                if ~ok_alejado180 || ~ok_alejado60 || ~ok_pen || ~ok_ult
                    fprintf('  RS->FA omitida t=%.2f s: alguna ventana no válida.\n', t_ini_FA);
                    continue
                end

                %% --------------------------------------------------------
                % Frecuencia dominante en cada ventana
                %% --------------------------------------------------------

                [DF_alejado180,~,~,Pow_alejado180] = frecuencia_dominante2_RS(seg_alejado180, Fs);
                [DF_alejado60,~,~,Pow_alejado60]   = frecuencia_dominante2_RS(seg_alejado60, Fs);
                [DF_pen,~,~,Pow_pen]               = frecuencia_dominante2_RS(seg_penult, Fs);
                [DF_ult,~,~,Pow_ult]               = frecuencia_dominante2_RS(seg_ult, Fs);

                if ~isfinite(DF_alejado180) || ~isfinite(DF_alejado60) || ...
                   ~isfinite(DF_pen) || ~isfinite(DF_ult)
                    continue
                end

                n_RS_FA_reg = n_RS_FA_reg + 1;

                resultados_RS_FA(end+1,:) = { ...
                    nombre_base, nombre_registro, ...
                    t_ini_FA, ...
                    DF_alejado180, DF_alejado60, DF_pen, DF_ult, ...
                    Pow_alejado180, Pow_alejado60, Pow_pen, Pow_ult, ...
                    numel(locs_qrs_local), ...
                    numel(locs_R_local)}; %#ok<SAGROW>

            end

            fprintf('  Transiciones RS->FA válidas: %d\n', n_RS_FA_reg);

        catch ME

            fprintf('  ERROR en %s: %s\n', nombre_registro, ME.message);

            if ~isempty(ME.stack)
                fprintf('  Línea aproximada: %d\n', ME.stack(1).line);
            end
        end
    end
end

%% ============================================================
% GUARDAR RESULTADOS RS -> FA
%% ============================================================

if ~isempty(resultados_RS_FA)

    T_RS_FA = cell2table(resultados_RS_FA, 'VariableNames', { ...
        'Base','Registro','Tiempo_ini_FA', ...
        'DF_alejado_m180_m175', ...
        'DF_alejado_m60_m55', ...
        'DF_penultimo_m10_m5', ...
        'DF_ultimo_m5_0', ...
        'Power_alejado_m180_m175', ...
        'Power_alejado_m60_m55', ...
        'Power_penultimo_m10_m5', ...
        'Power_ultimo_m5_0', ...
        'N_R_qrs_local', ...
        'N_R_hibrido_local'});

    T_RS_FA = T_RS_FA( ...
        isfinite(T_RS_FA.DF_alejado_m180_m175) & ...
        isfinite(T_RS_FA.DF_alejado_m60_m55) & ...
        isfinite(T_RS_FA.DF_penultimo_m10_m5) & ...
        isfinite(T_RS_FA.DF_ultimo_m5_0), :);

    if ~isempty(T_RS_FA)

        %% Deltas

        T_RS_FA.Delta_DF_ultimo_menos_penultimo = ...
            T_RS_FA.DF_ultimo_m5_0 - T_RS_FA.DF_penultimo_m10_m5;

        T_RS_FA.Delta_DF_ultimo_menos_60s = ...
            T_RS_FA.DF_ultimo_m5_0 - T_RS_FA.DF_alejado_m60_m55;

        T_RS_FA.Delta_DF_ultimo_menos_180s = ...
            T_RS_FA.DF_ultimo_m5_0 - T_RS_FA.DF_alejado_m180_m175;

        writetable(T_RS_FA, fullfile(carpeta_out,'resultados_RS_FA_ventanas.xlsx'));

        %% Estadísticos descriptivos

        media_180 = mean(T_RS_FA.DF_alejado_m180_m175, 'omitnan');
        std_180   = std(T_RS_FA.DF_alejado_m180_m175, 0, 'omitnan');

        media_60 = mean(T_RS_FA.DF_alejado_m60_m55, 'omitnan');
        std_60   = std(T_RS_FA.DF_alejado_m60_m55, 0, 'omitnan');

        media_pen = mean(T_RS_FA.DF_penultimo_m10_m5, 'omitnan');
        std_pen   = std(T_RS_FA.DF_penultimo_m10_m5, 0, 'omitnan');

        media_ult = mean(T_RS_FA.DF_ultimo_m5_0, 'omitnan');
        std_ult   = std(T_RS_FA.DF_ultimo_m5_0, 0, 'omitnan');

        media_delta_pen = mean(T_RS_FA.Delta_DF_ultimo_menos_penultimo, 'omitnan');
        std_delta_pen   = std(T_RS_FA.Delta_DF_ultimo_menos_penultimo, 0, 'omitnan');

        media_delta_60 = mean(T_RS_FA.Delta_DF_ultimo_menos_60s, 'omitnan');
        std_delta_60   = std(T_RS_FA.Delta_DF_ultimo_menos_60s, 0, 'omitnan');

        media_delta_180 = mean(T_RS_FA.Delta_DF_ultimo_menos_180s, 'omitnan');
        std_delta_180   = std(T_RS_FA.Delta_DF_ultimo_menos_180s, 0, 'omitnan');

        %% Wilcoxon pareado

        if height(T_RS_FA) >= 2

            p_wilcoxon_ultimo_vs_penultimo = signrank( ...
                T_RS_FA.DF_ultimo_m5_0, ...
                T_RS_FA.DF_penultimo_m10_m5);

            p_wilcoxon_ultimo_vs_60s = signrank( ...
                T_RS_FA.DF_ultimo_m5_0, ...
                T_RS_FA.DF_alejado_m60_m55);

            p_wilcoxon_ultimo_vs_180s = signrank( ...
                T_RS_FA.DF_ultimo_m5_0, ...
                T_RS_FA.DF_alejado_m180_m175);

        else
            p_wilcoxon_ultimo_vs_penultimo = NaN;
            p_wilcoxon_ultimo_vs_60s = NaN;
            p_wilcoxon_ultimo_vs_180s = NaN;
        end

        tabla_resumen_RS_FA = table( ...
            height(T_RS_FA), ...
            {sprintf('%.3f +/- %.3f', media_180, std_180)}, ...
            {sprintf('%.3f +/- %.3f', media_60, std_60)}, ...
            {sprintf('%.3f +/- %.3f', media_pen, std_pen)}, ...
            {sprintf('%.3f +/- %.3f', media_ult, std_ult)}, ...
            {sprintf('%.3f +/- %.3f', media_delta_pen, std_delta_pen)}, ...
            {sprintf('%.3f +/- %.3f', media_delta_60, std_delta_60)}, ...
            {sprintf('%.3f +/- %.3f', media_delta_180, std_delta_180)}, ...
            {sprintf('%.2e', p_wilcoxon_ultimo_vs_penultimo)}, ...
            {sprintf('%.2e', p_wilcoxon_ultimo_vs_60s)}, ...
            {sprintf('%.2e', p_wilcoxon_ultimo_vs_180s)}, ...
            'VariableNames', { ...
            'N_transiciones', ...
            'DF_alejado_m180_m175_Hz', ...
            'DF_alejado_m60_m55_Hz', ...
            'DF_penultimo_m10_m5_Hz', ...
            'DF_ultimo_m5_0_Hz', ...
            'Delta_ultimo_menos_penultimo_Hz', ...
            'Delta_ultimo_menos_60s_Hz', ...
            'Delta_ultimo_menos_180s_Hz', ...
            'p_wilcoxon_ultimo_vs_penultimo', ...
            'p_wilcoxon_ultimo_vs_60s', ...
            'p_wilcoxon_ultimo_vs_180s'});

        writetable(tabla_resumen_RS_FA, ...
            fullfile(carpeta_out,'resumen_estadistico_RS_FA_ventanas.xlsx'));

        disp(' ');
        disp('RESUMEN ESTADISTICO RS -> FA');
        disp(tabla_resumen_RS_FA);

        %% ============================================================
        % FIGURAS
        %% ============================================================

        %% Boxplot de las cuatro ventanas

        f1 = figure('Visible','off','Color','w', 'Position', [100 100 1100 650]);

        boxplot([ ...
            T_RS_FA.DF_alejado_m180_m175, ...
            T_RS_FA.DF_alejado_m60_m55, ...
            T_RS_FA.DF_penultimo_m10_m5, ...
            T_RS_FA.DF_ultimo_m5_0], ...
            {'-180 a -175 s','-60 a -55 s','-10 a -5 s','-5 a 0 s'});

        title('Frecuencia dominante del residual auricular previa al inicio de la FA', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f1, fullfile(carpeta_out, 'boxplot_DF_RS_FA_4ventanas.png'));
        close(f1);

        %% Boxplot cercano: -10:-5 vs -5:0

        f2 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        boxplot([ ...
            T_RS_FA.DF_penultimo_m10_m5, ...
            T_RS_FA.DF_ultimo_m5_0], ...
            {'-10 a -5 s','-5 a 0 s'});

        title('DF residual auricular: comparación final antes del inicio de la FA', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f2, fullfile(carpeta_out, 'boxplot_DF_RS_FA_penultimo_vs_ultimo.png'));
        close(f2);

        %% Boxplot 60 s: -60:-55 vs -5:0

        f3 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        boxplot([ ...
            T_RS_FA.DF_alejado_m60_m55, ...
            T_RS_FA.DF_ultimo_m5_0], ...
            {'-60 a -55 s','-5 a 0 s'});

        title('DF residual auricular: ventana -60 s frente a ventana final', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f3, fullfile(carpeta_out, 'boxplot_DF_RS_FA_60s_vs_ultimo.png'));
        close(f3);

        %% Boxplot 180 s: -180:-175 vs -5:0

        f4 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        boxplot([ ...
            T_RS_FA.DF_alejado_m180_m175, ...
            T_RS_FA.DF_ultimo_m5_0], ...
            {'-180 a -175 s','-5 a 0 s'});

        title('DF residual auricular: ventana -180 s frente a ventana final', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f4, fullfile(carpeta_out, 'boxplot_DF_RS_FA_180s_vs_ultimo.png'));
        close(f4);

        %% Gráfico pareado: ventana -180:-175 frente a última ventana

        f_pareado2 = figure('Visible','off','Color','w', 'Position', [100 100 1100 650]);
        
        DF_alej = T_RS_FA.DF_alejado_m180_m175;
        DF_ult  = T_RS_FA.DF_ultimo_m5_0;
        
        idx_validos = isfinite(DF_alej) & isfinite(DF_ult);
        
        DF_alej = DF_alej(idx_validos);
        DF_ult  = DF_ult(idx_validos);
        
        hold on
        
        for k = 1:numel(DF_alej)

            plot([1 2], [DF_alej(k), DF_ult(k)], ...
                'Color', [0.80 0.80 0.80], ...
                'LineWidth', 0.6);
        
            plot(1, DF_alej(k), ...
                'o', ...
                'Color', [0.35 0.35 0.35], ...
                'MarkerSize', 6);
        
            plot(2, DF_ult(k), ...
                'o', ...
                'Color', [0.35 0.35 0.35], ...
                'MarkerSize', 6);
        end

        % Línea de la media
        plot([1 2], [mean(DF_alej), mean(DF_ult)], ...
            'k-', ...
            'LineWidth', 2.2);

        plot(1, mean(DF_alej), ...
            'ks', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 8);

        plot(2, mean(DF_ult), ...
            'ks', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 8);
        
        hold off
        
        xlim([0.6 2.4])
        ylim([0.5 4])
        
        set(gca, ...
            'XTick', [1 2], ...
            'XTickLabel', {'-180 a -175 s', '-5 a 0 s'}, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        ylabel('Frecuencia dominante (Hz)', ...
            'Interpreter','none', ...
            'FontSize', 20);

        xlabel('Ventana temporal', ...
            'Interpreter','none', ...
            'FontSize', 20);

        title('Evolución pareada de la DF: -180 a -175 s frente a -5 a 0 s', ...
            'Interpreter','none', ...
            'FontSize', 18, ...
            'FontWeight','bold');
        
        grid on
        box on
        
        saveas(f_pareado2, fullfile(carpeta_out, ...
            'grafico_pareado_DF_RS_FA_alejado180_vs_ultimo.png'));
        
        close(f_pareado2);

        %% Gráfico pareado: penúltima frente a última ventana

        f_pareado1 = figure('Visible','off','Color','w', 'Position', [100 100 1100 650]);

        DF_pen = T_RS_FA.DF_penultimo_m10_m5;
        DF_ult = T_RS_FA.DF_ultimo_m5_0;

        idx_validos = isfinite(DF_pen) & isfinite(DF_ult);

        DF_pen = DF_pen(idx_validos);
        DF_ult = DF_ult(idx_validos);

        hold on

        for k = 1:numel(DF_pen)

            plot([1 2], [DF_pen(k), DF_ult(k)], ...
                'Color', [0.80 0.80 0.80], ...
                'LineWidth', 0.6);

            plot(1, DF_pen(k), ...
                'o', ...
                'Color', [0.35 0.35 0.35], ...
                'MarkerSize', 6);

            plot(2, DF_ult(k), ...
                'o', ...
                'Color', [0.35 0.35 0.35], ...
                'MarkerSize', 6);
        end

        % Línea de la media
        plot([1 2], [mean(DF_pen), mean(DF_ult)], ...
            'k-', ...
            'LineWidth', 2.2);

        plot(1, mean(DF_pen), ...
            'ks', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 8);

        plot(2, mean(DF_ult), ...
            'ks', ...
            'MarkerFaceColor', 'k', ...
            'MarkerSize', 8);

        hold off

        xlim([0.6 2.4])
        ylim([0.5 4])

        set(gca, ...
            'XTick', [1 2], ...
            'XTickLabel', {'-10 a -5 s', '-5 a 0 s'}, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        ylabel('Frecuencia dominante (Hz)', ...
            'Interpreter','none', ...
            'FontSize', 20);

        xlabel('Ventana temporal', ...
            'Interpreter','none', ...
            'FontSize', 20);

        title('Evolución pareada de la DF: -10 a -5 frente a -5 a 0', ...
            'Interpreter','none', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        grid on
        box on

        saveas(f_pareado1, fullfile(carpeta_out, ...
            'grafico_pareado_DF_RS_FA_penultimo_vs_ultimo.png'));

        close(f_pareado1);
        
        %% Histograma delta cercano

        f5 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        histogram(T_RS_FA.Delta_DF_ultimo_menos_penultimo);
        xline(0, '--', 'LineWidth', 1.2);

        title('Cambio de DF: (-5 a 0 s) - (-10 a -5 s)', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        xlabel('\Delta DF (Hz)', ...
            'FontSize', 20);

        ylabel('Número de transiciones', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f5, fullfile(carpeta_out, 'histograma_Delta_DF_ultimo_menos_penultimo.png'));
        close(f5);

        %% Histograma delta 60 s

        f6 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        histogram(T_RS_FA.Delta_DF_ultimo_menos_60s);
        xline(0, '--', 'LineWidth', 1.2);

        title('Cambio de DF: (-5 a 0 s) - (-60 a -55 s)', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        xlabel('\Delta DF (Hz)', ...
            'FontSize', 20);

        ylabel('Número de transiciones', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f6, fullfile(carpeta_out, 'histograma_Delta_DF_ultimo_menos_60s.png'));
        close(f6);

        %% Histograma delta 180 s

        f7 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        histogram(T_RS_FA.Delta_DF_ultimo_menos_180s);
        xline(0, '--', 'LineWidth', 1.2);

        title('Cambio de DF: (-5 a 0 s) - (-180 a -175 s)', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        xlabel('\Delta DF (Hz)', ...
            'FontSize', 20);

        ylabel('Número de transiciones', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f7, fullfile(carpeta_out, 'histograma_Delta_DF_ultimo_menos_180s.png'));
        close(f7);

    else
        warning('Resultados RS->FA existen pero todos son NaN tras limpiar.');
    end

else
    warning('No se encontraron resultados RS -> FA.');
end

fprintf('\nFIN\n');
fprintf('Resultados guardados en:\n%s\n', carpeta_out);

%% ============================================================
% FUNCIONES LOCALES
%% ============================================================

function [segmento, ok] = extraer_ventana_relativa_fin(x, Fs, t_ini_rel, t_fin_rel)

% Extrae una ventana definida respecto al final de x.
%
% En este script, x corresponde a los 180 s previos al inicio de la FA,
% por lo que el final de x coincide con t_ini_FA.
%
% Ejemplo:
%   t_ini_rel = -180
%   t_fin_rel = -175
%
% Devuelve el segmento situado entre 180 y 175 s antes del inicio de la FA.

segmento = [];
ok = false;

x = x(:);
N = length(x);

if isempty(x) || Fs <= 0 || t_fin_rel <= t_ini_rel
    return
end

dur_total = N / Fs;

t_ini_abs = dur_total + t_ini_rel;
t_fin_abs = dur_total + t_fin_rel;

idx_ini = round(t_ini_abs * Fs) + 1;
idx_fin = round(t_fin_abs * Fs);

if idx_ini < 1 || idx_fin > N || idx_fin <= idx_ini
    return
end

segmento = x(idx_ini:idx_fin);
segmento = segmento(:);

dur_esperada = round((t_fin_rel - t_ini_rel) * Fs);

if length(segmento) ~= dur_esperada
    segmento = segmento(1:min(end,dur_esperada));
end

if length(segmento) < dur_esperada
    return
end

ok = true;

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

function ritmo = mapear_ritmo(txt)
    if strcmp(txt,'(AFIB')
        ritmo = 'FA';
    elseif strcmp(txt,'(N') || strcmp(txt,'(NSR') || strcmp(txt,'(SR')
        ritmo = 'RS';
    else
        ritmo = 'OTHER';
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

function [x_residual, plantilla, t_plantilla, latidos_validos] = cancelar_QRST_plantilla_medianaRS(x, locs_R, Fs)

% CANCELAR_QRST_PLANTILLA_MEDIANARS
% Cancela los complejos QRS-T mediante sustracción de una plantilla mediana.
%
% Versión para RS:
%   - 60 ms antes del R
%   - 450 ms después del R
%
% Entrada:
%   x      -> señal ECG filtrada
%   locs_R -> posiciones de los picos R en muestras
%   Fs     -> frecuencia de muestreo
%
% Salida:
%   x_residual      -> señal residual con menor contribución QRS-T
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
post_R = round(0.45 * Fs);   % 450 ms después del R

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

    % Ajuste de amplitud para adaptar la plantilla a cada latido
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

%% 5) CENTRAR SEÑAL RESIDUAL

x_residual = x_residual - mean(x_residual, 'omitnan');

end

function [DF, f_axis, Pxx, peak_power] = frecuencia_dominante2_RS(x, Fs)

DF = NaN;
peak_power = NaN;
f_axis = [];
Pxx = [];

x = x(:);

if isempty(x) || numel(x) < Fs
    return
end

if any(~isfinite(x))
    return
end

%% Quitar media y tendencia

x = x - mean(x, 'omitnan');
x = detrend(x);

%% Espectro de potencia con pwelch

N = length(x);
Nfft = 8192;

ventana = hamming(N);
noverlap = 0;

[Pxx, f_axis] = pwelch(x, ventana, noverlap, Nfft, Fs);

%% Buscar pico dominante entre 0.5 y 2 Hz

idx = f_axis >= 0.5 & f_axis <= 2;

if ~any(idx)
    return
end

P_banda = Pxx(idx);
f_valid = f_axis(idx);

[Pmax, im] = max(P_banda);

DF = f_valid(im);
peak_power = P_banda(im);

%% Control de calidad: pico dominante claro

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

%% Rechazar picos pegados a los bordes

if DF <= 0.55 || DF >= 1.95
    DF = NaN;
    peak_power = NaN;
    return
end

end
