clear
clc
close all

%% ============================================================
% ANÁLISIS DE TRANSICIONES FA -> RS
%
% Este script analiza transiciones de fibrilación auricular a ritmo sinusal
% en registros de larga duración. Para cada transición válida se estudia la
% frecuencia dominante (DF) del residual auricular en tres ventanas de 2 s
% previas al final del episodio de FA.
%
% Estrategia general:
%   - Se leen la señal ECG y las anotaciones de ritmo de los registros.
%   - Se descarta el registro 112.
%   - Se detectan y excluyen ventanas con artefactos.
%   - Se utilizan las anotaciones .qrs como referencia principal para los
%     picos R.
%   - DetectarPICOSR complementa los picos R cuando existen huecos largos.
%   - Se realiza un ajuste final de los picos R al máximo local.
%   - Se cancela el complejo QRS-T mediante plantilla mediana.
%   - Se calcula la DF del residual auricular en ventanas próximas al fin
%     de la FA.
%
% Ventanas analizadas respecto al final de la FA:
%   - Ventana alejada:   -12 a -10 s
%   - Ventana penúltima:  -4 a  -2 s
%   - Ventana última:     -2 a   0 s
%
% Comparaciones realizadas:
%   - Delta clásico = DF(-2 a 0 s) - DF(-4 a -2 s)
%   - Delta alejado = DF(-2 a 0 s) - DF(-12 a -10 s)
%
% Salidas:
%   - resultados_FA_RS_ventanas.xlsx
%   - resumen_estadistico_FA_RS_ventanas.xlsx
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
%   - cancelar_QRST_plantilla_medianaFA
%   - frecuencia_dominante2_FA
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

carpeta_out = 'C:\Users\Emma\Documents\MATLAB\Petrutiu_FA_RS_sinartefactos_definitivo';

if ~exist(carpeta_out, 'dir')
    mkdir(carpeta_out);
end

Fs_new = 500;

dur_cancel_FA = 30;   % segundos antes del fin de FA para cancelar QRST

vent_alejada_ini = -12;
vent_alejada_fin = -10;

vent_pen_ini = -4;
vent_pen_fin = -2;

vent_ult_ini = -2;
vent_ult_fin = 0;

min_R_FA_cancel = 15;

registros_excluir = {'112'};

resultados_FA_RS = {};

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
            fprintf('\nRegistro %s excluido del analisis.\n', nombre_registro);
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
                fprintf('  Senal vacia. Se omite.\n');
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
            % El filtrado 1-50 Hz del análisis de FA se mantiene sin cambios.
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
            % 0 buena, 1 desconexion, 2 transicion, 3 artefacto

            %% Detectar desconexion

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

            %% Detectar transicion despues de desconexion

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
            % 3) FILTRADO PARA ANÁLISIS DE FA
            %% ============================================================

            [b_FA, a_FA] = butter(2, [1 50] / (Fs/2), 'bandpass');
            x_filt_FA = filtfilt(b_FA, a_FA, x);

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
            locs_R_qrs_global = locs_R_qrs_global(locs_R_qrs_global >= 1 & ...
                                                  locs_R_qrs_global <= length(x));

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
            % Las anotaciones consecutivas que representen el mismo ritmo
            % se unifican despues de ordenarlas, conservando la primera
            % anotacion de cada episodio continuo.
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
                fprintf('  No hay ritmos validos. Se omite.\n');
                continue
            end

            [tiempos_ritmo, idx] = sort(tiempos_ritmo);
            ritmos = ritmos(idx);

            %% ============================================================
            % 6.1) UNIFICAR ANOTACIONES CONSECUTIVAS DEL MISMO RITMO
            %
            % Si aparecen varias anotaciones consecutivas que representan
            % el mismo ritmo, se consideran parte de un unico episodio.
            %
            % Ejemplo:
            %   0 s    (AFIB
            %   400 s  (AFIB
            %   500 s  (N
            %
            % Se transforma en:
            %   0 s    (AFIB
            %   500 s  (N
            %
            % Asi se conserva el inicio real del episodio continuo de FA.

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
                    'tras unificar repeticiones.\n'])
                continue
            end

            %% ============================================================
            % 7) ANÁLISIS FA -> RS
            %% ============================================================

            n_FA_RS_reg = 0;

            for i = 1:length(tiempos_ritmo)-1

                ritmo_actual = mapear_ritmo(ritmos{i});
                ritmo_sig    = mapear_ritmo(ritmos{i+1});

                % Compatibilidad por si mapear_ritmo devuelve AF/SR o FA/RS
                es_FA_actual = strcmp(ritmo_actual,'FA') || strcmp(ritmo_actual,'AF');
                es_RS_sig    = strcmp(ritmo_sig,'RS')    || strcmp(ritmo_sig,'SR');

                if ~(es_FA_actual && es_RS_sig)
                    continue
                end

                t_ini_FA = tiempos_ritmo(i);
                t_fin_FA = tiempos_ritmo(i+1);

                if (t_fin_FA - t_ini_FA) < dur_cancel_FA
                    continue
                end

                %% --------------------------------------------------------
                % Comprobar artefactos en las tres ventanas analizadas
                %
                % Solo se descarta la transicion si un intervalo malo se
                % solapa con:
                %   -12 a -10 s
                %    -4 a  -2 s
                %    -2 a   0 s
                %% --------------------------------------------------------

                ventanas_calidad = [ ...
                    t_fin_FA + vent_alejada_ini, ...
                    t_fin_FA + vent_alejada_fin; ...
                    t_fin_FA + vent_pen_ini, ...
                    t_fin_FA + vent_pen_fin; ...
                    t_fin_FA + vent_ult_ini, ...
                    t_fin_FA + vent_ult_fin];

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

                    fprintf(['  FA->RS omitida t=%.2f s: alguna ventana ' ...
                        'analizada contiene artefactos.\n'], t_fin_FA);

                    continue
                end

                %% --------------------------------------------------------
                % Ventana local de 30 s antes del fin de FA
                %% --------------------------------------------------------

                t0_local = t_fin_FA - dur_cancel_FA;
                t1_local = t_fin_FA;

                ini_local = round(t0_local * Fs) + 1;
                fin_local = round(t1_local * Fs);

                if ini_local < 1 || fin_local > length(x_filt_FA)
                    continue
                end

                ecg_local_FA = x_filt_FA(ini_local:fin_local);
                ecg_local_FA = ecg_local_FA(:);

                if length(ecg_local_FA) < dur_cancel_FA * Fs
                    continue
                end

                %% --------------------------------------------------------
                % R .qrs locales
                %% --------------------------------------------------------

                idx_qrs_local = locs_R_qrs_global >= ini_local & locs_R_qrs_global <= fin_local;
                locs_qrs_local = locs_R_qrs_global(idx_qrs_local) - ini_local + 1;
                locs_qrs_local = limpiar_locs_local(locs_qrs_local, length(ecg_local_FA));

                %% --------------------------------------------------------
                % R hibridos locales
                %% --------------------------------------------------------

                [locs_R_local, motivo_R] = DetectarPICOSR(ecg_local_FA, locs_qrs_local, Fs);

                if isempty(locs_R_local) || numel(locs_R_local) < min_R_FA_cancel
                    fprintf('  FA->RS omitida t=%.2f s: R insuficientes. %s\n', ...
                        t_fin_FA, motivo_R);
                    continue
                end

                %% --------------------------------------------------------
                % Cancelar QRST localmente
                %% --------------------------------------------------------

                atrial_local_FA = cancelar_QRST_plantilla_medianaFA(ecg_local_FA, locs_R_local, Fs);

                if isempty(atrial_local_FA) || any(~isfinite(atrial_local_FA))
                    fprintf('  FA->RS omitida t=%.2f s: residual FA no valido.\n', t_fin_FA);
                    continue
                end

                atrial_local_FA = atrial_local_FA(:);

                %% --------------------------------------------------------
                % Extraer las tres ventanas respecto al final de FA
                %
                % Como atrial_local_FA va desde:
                %   t_fin_FA - 30 s  hasta  t_fin_FA
                %
                % El instante relativo 0 corresponde al final de la FA.
                %% --------------------------------------------------------

                [seg_alejado, ok_alejado] = extraer_ventana_relativa_fin( ...
                    atrial_local_FA, Fs, vent_alejada_ini, vent_alejada_fin);

                [seg_penult, ok_pen] = extraer_ventana_relativa_fin( ...
                    atrial_local_FA, Fs, vent_pen_ini, vent_pen_fin);

                [seg_ult, ok_ult] = extraer_ventana_relativa_fin( ...
                    atrial_local_FA, Fs, vent_ult_ini, vent_ult_fin);

                if ~ok_alejado || ~ok_pen || ~ok_ult
                    fprintf('  FA->RS omitida t=%.2f s: alguna ventana no valida.\n', t_fin_FA);
                    continue
                end

                %% --------------------------------------------------------
                % Frecuencia dominante en cada ventana
                %% --------------------------------------------------------

                [DF_alejado,~,~,Pow_alejado] = frecuencia_dominante2_FA(seg_alejado, Fs);
                [DF_pen,~,~,Pow_pen]         = frecuencia_dominante2_FA(seg_penult, Fs);
                [DF_ult,~,~,Pow_ult]         = frecuencia_dominante2_FA(seg_ult, Fs);

                if ~isfinite(DF_alejado) || ~isfinite(DF_pen) || ~isfinite(DF_ult)
                    continue
                end

                n_FA_RS_reg = n_FA_RS_reg + 1;

                resultados_FA_RS(end+1,:) = { ...
                    nombre_base, nombre_registro, ...
                    t_ini_FA, t_fin_FA, ...
                    DF_alejado, DF_pen, DF_ult, ...
                    Pow_alejado, Pow_pen, Pow_ult, ...
                    numel(locs_qrs_local), ...
                    numel(locs_R_local)}; %#ok<SAGROW>

            end

            fprintf('  Transiciones FA->RS validas: %d\n', n_FA_RS_reg);

        catch ME

            fprintf('  ERROR en %s: %s\n', nombre_registro, ME.message);

            if ~isempty(ME.stack)
                fprintf('  Linea aproximada: %d\n', ME.stack(1).line);
            end
        end
    end
end

%% ============================================================
% GUARDAR RESULTADOS FA -> RS
%% ============================================================

if ~isempty(resultados_FA_RS)

    T_FA_RS = cell2table(resultados_FA_RS, 'VariableNames', { ...
        'Base','Registro','Tiempo_ini_FA','Tiempo_fin_FA', ...
        'DF_alejado_m12_m10', ...
        'DF_penultimo_m4_m2', ...
        'DF_ultimo_m2_0', ...
        'Power_alejado_m12_m10', ...
        'Power_penultimo_m4_m2', ...
        'Power_ultimo_m2_0', ...
        'N_R_qrs_local', ...
        'N_R_hibrido_local'});

    T_FA_RS = T_FA_RS( ...
        isfinite(T_FA_RS.DF_alejado_m12_m10) & ...
        isfinite(T_FA_RS.DF_penultimo_m4_m2) & ...
        isfinite(T_FA_RS.DF_ultimo_m2_0), :);

    if ~isempty(T_FA_RS)

        %% Deltas

        T_FA_RS.Delta_DF_ultimo_menos_penultimo = ...
            T_FA_RS.DF_ultimo_m2_0 - T_FA_RS.DF_penultimo_m4_m2;

        T_FA_RS.Delta_DF_ultimo_menos_alejado = ...
            T_FA_RS.DF_ultimo_m2_0 - T_FA_RS.DF_alejado_m12_m10;

        writetable(T_FA_RS, fullfile(carpeta_out,'resultados_FA_RS_ventanas.xlsx'));

        %% Estadísticos descriptivos

        media_alejado = mean(T_FA_RS.DF_alejado_m12_m10, 'omitnan');
        std_alejado   = std(T_FA_RS.DF_alejado_m12_m10, 0, 'omitnan');

        media_pen = mean(T_FA_RS.DF_penultimo_m4_m2, 'omitnan');
        std_pen   = std(T_FA_RS.DF_penultimo_m4_m2, 0, 'omitnan');

        media_ult = mean(T_FA_RS.DF_ultimo_m2_0, 'omitnan');
        std_ult   = std(T_FA_RS.DF_ultimo_m2_0, 0, 'omitnan');

        media_delta_pen = mean(T_FA_RS.Delta_DF_ultimo_menos_penultimo, 'omitnan');
        std_delta_pen   = std(T_FA_RS.Delta_DF_ultimo_menos_penultimo, 0, 'omitnan');

        media_delta_alejado = mean(T_FA_RS.Delta_DF_ultimo_menos_alejado, 'omitnan');
        std_delta_alejado   = std(T_FA_RS.Delta_DF_ultimo_menos_alejado, 0, 'omitnan');

        %% Wilcoxon pareado

        if height(T_FA_RS) >= 2

            p_wilcoxon_ultimo_vs_penultimo = signrank( ...
                T_FA_RS.DF_ultimo_m2_0, ...
                T_FA_RS.DF_penultimo_m4_m2);

            p_wilcoxon_ultimo_vs_alejado = signrank( ...
                T_FA_RS.DF_ultimo_m2_0, ...
                T_FA_RS.DF_alejado_m12_m10);

        else

            p_wilcoxon_ultimo_vs_penultimo = NaN;
            p_wilcoxon_ultimo_vs_alejado = NaN;

        end

        tabla_resumen_FA_RS = table( ...
            height(T_FA_RS), ...
            {sprintf('%.3f +/- %.3f', media_alejado, std_alejado)}, ...
            {sprintf('%.3f +/- %.3f', media_pen, std_pen)}, ...
            {sprintf('%.3f +/- %.3f', media_ult, std_ult)}, ...
            {sprintf('%.3f +/- %.3f', media_delta_pen, std_delta_pen)}, ...
            {sprintf('%.3f +/- %.3f', media_delta_alejado, std_delta_alejado)}, ...
            {sprintf('%.2e', p_wilcoxon_ultimo_vs_penultimo)}, ...
            {sprintf('%.2e', p_wilcoxon_ultimo_vs_alejado)}, ...
            'VariableNames', { ...
            'N_transiciones', ...
            'DF_alejado_m12_m10_Hz', ...
            'DF_penultimo_m4_m2_Hz', ...
            'DF_ultimo_m2_0_Hz', ...
            'Delta_ultimo_menos_penultimo_Hz', ...
            'Delta_ultimo_menos_alejado_Hz', ...
            'p_wilcoxon_ultimo_vs_penultimo', ...
            'p_wilcoxon_ultimo_vs_alejado'});

        writetable(tabla_resumen_FA_RS, ...
            fullfile(carpeta_out,'resumen_estadistico_FA_RS_ventanas.xlsx'));

        disp(' ');
        disp('RESUMEN ESTADISTICO FA -> RS');
        disp(tabla_resumen_FA_RS);

        %% ============================================================
        % FIGURAS
        %% ============================================================

        %% Boxplot de las tres ventanas

        f1 = figure('Visible','off','Color','w', 'Position', [100 100 1100 650]);

        boxplot([ ...
            T_FA_RS.DF_alejado_m12_m10, ...
            T_FA_RS.DF_penultimo_m4_m2, ...
            T_FA_RS.DF_ultimo_m2_0], ...
            {'-12 a -10 s','-4 a -2 s','-2 a 0 s'});

        title('Frecuencia dominante fibrilatoria previa a la terminacion de la FA', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f1, fullfile(carpeta_out, 'boxplot_DF_FA_RS_3ventanas.png'));
        close(f1);

        %% Boxplot clasico: -4:-2 vs -2:0

        f2 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        boxplot([ ...
            T_FA_RS.DF_penultimo_m4_m2, ...
            T_FA_RS.DF_ultimo_m2_0], ...
            {'-4 a -2 s','-2 a 0 s'});

        title('DF fibrilatoria: comparacion final antes de la terminacion de la FA', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f2, fullfile(carpeta_out, 'boxplot_DF_FA_RS_penultimo_vs_ultimo.png'));
        close(f2);

        %% Boxplot alejado: -12:-10 vs -2:0

        f3 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        boxplot([ ...
            T_FA_RS.DF_alejado_m12_m10, ...
            T_FA_RS.DF_ultimo_m2_0], ...
            {'-12 a -10 s','-2 a 0 s'});

        title('DF fibrilatoria: ventana alejada frente a ventana final', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f3, fullfile(carpeta_out, 'boxplot_DF_FA_RS_alejado_vs_ultimo.png'));
        close(f3);

        %% Gráfico pareado: ventana alejada frente a última ventana

        f_pareado1 = figure('Visible','off','Color','w', 'Position', [100 100 1100 650]);

        DF_alej = T_FA_RS.DF_alejado_m12_m10;
        DF_ult  = T_FA_RS.DF_ultimo_m2_0;

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

        set(gca, ...
            'XTick', [1 2], ...
            'XTickLabel', {'-12 a -10 s','-2 a 0 s'}, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        xlabel('Ventana temporal', ...
            'FontSize', 20);

        title('Evolucion pareada de la DF: -12 a -10 s frente a -2 a 0 s ', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        grid on
        box on

        saveas(f_pareado1, fullfile(carpeta_out, ...
            'grafico_pareado_DF_FA_RS_alejado_vs_ultimo.png'));

        close(f_pareado1);

        %% Gráfico pareado clásico: penúltima vs última ventana

        f_pareado2 = figure('Visible','off','Color','w', 'Position', [100 100 1100 650]);

        DF_pen = T_FA_RS.DF_penultimo_m4_m2;
        DF_ult = T_FA_RS.DF_ultimo_m2_0;

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

        set(gca, ...
            'XTick', [1 2], ...
            'XTickLabel', {'-4 a -2 s','-2 a 0 s'}, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        ylabel('Frecuencia dominante (Hz)', ...
            'FontSize', 20);

        xlabel('Ventana temporal', ...
            'FontSize', 20);

        title('Evolucion pareada de la DF: -4 a -2 s frente a -2 a 0 s', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        grid on
        box on

        saveas(f_pareado2, fullfile(carpeta_out, ...
            'grafico_pareado_DF_FA_RS_penultimo_vs_ultimo.png'));

        close(f_pareado2);

        %% Histograma delta clásico

        f4 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        histogram(T_FA_RS.Delta_DF_ultimo_menos_penultimo);

        title('Cambio de DF: (-2 a 0 s) - (-4 a -2 s)', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        xlabel('\Delta DF (Hz)', ...
            'FontSize', 20);

        ylabel('Numero de transiciones', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f4, fullfile(carpeta_out, ...
            'histograma_Delta_DF_ultimo_menos_penultimo.png'));

        close(f4);

        %% Histograma delta alejado

        f5 = figure('Visible','off','Color','w', 'Position', [100 100 950 650]);

        histogram(T_FA_RS.Delta_DF_ultimo_menos_alejado);

        title('Cambio de DF: (-2 a 0 s) - (-12 a -10 s)', ...
            'FontSize', 18, ...
            'FontWeight','bold');

        xlabel('\Delta DF (Hz)', ...
            'FontSize', 20);

        ylabel('Numero de transiciones', ...
            'FontSize', 20);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas(f5, fullfile(carpeta_out, ...
            'histograma_Delta_DF_ultimo_menos_alejado.png'));

        close(f5);

    else
        warning('Resultados FA->RS existen pero todos son NaN tras limpiar.');
    end

else
    warning('No se encontraron resultados FA -> RS.');
end

fprintf('\nFIN\n');
fprintf('Resultados guardados en:\n%s\n', carpeta_out);

%% ============================================================
% FUNCIONES LOCALES
%% ============================================================

function [segmento, ok] = extraer_ventana_relativa_fin(x, Fs, t_ini_rel, t_fin_rel)

% Extrae una ventana definida respecto al final de x.
%
% Ejemplo:
%   t_ini_rel = -12
%   t_fin_rel = -10
%
% Devuelve el segmento situado entre -12 y -10 s antes del final.

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

function [x_residual, plantilla, t_plantilla, latidos_validos] = cancelar_QRST_plantilla_medianaFA(x, locs_R, Fs)

% CANCELAR_QRST_PLANTILLA_MEDIANAFA
% Cancela los complejos QRS-T mediante sustracción de una plantilla mediana.
%
% Versión para FA:
%   - 200 ms antes del R
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
%   t_plantilla     -> eje temporal de la plantilla respecto al R, en segundos
%   latidos_validos -> latidos usados para construir la plantilla

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

function [DF, f_axis, Pxx, peak_power] = frecuencia_dominante2_FA(x, Fs)

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

%% Espectro de potencia con pwelch usando toda la señal
%
% No se vuelve a filtrar en 3-9 Hz.
% La señal de entrada ya debe ser el residual tras QRST,
% obtenido a partir del ECG filtrado correspondiente.

N = length(x);
Nfft = 8192;

ventana = hamming(N);
noverlap = 0;

[Pxx, f_axis] = pwelch(x, ventana, noverlap, Nfft, Fs);

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

if DF <= 3.1 || DF >= 8.9
    DF = NaN;
    peak_power = NaN;
    return
end

end