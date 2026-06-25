clear
clc
close all

%% ============================================================
% EXTRACCION DE BIOMARCADORES DE 10 s: TRANSICIONES RS -> FA
%
% Se leen directamente LTAFDB y AFDB.
%
% Para cada transicion RS -> FA se estudian los 180 s ANTERIORES,
% que deben pertenecer de forma continua a RS.
%
% La cancelacion QRST se realiza una sola vez sobre el bloque completo
% de 180 s. Despues se extraen del residual las tres ventanas de 10 s.
% Los biomarcadores RR y las ondas P/T se calculan de forma independiente
% dentro de cada ventana de 10 s.
%
% El control de artefactos se aplica unicamente a las tres ventanas
% de 10 s utilizadas para calcular los biomarcadores:
%
%   MOMENTO 1: -180 a -170 s  (RS alejado del inicio de FA)
%   MOMENTO 2:  -90 a  -80 s  (RS intermedio antes del inicio de FA)
%   MOMENTO 3:  -10 a    0 s  (RS inmediatamente previo a FA)
%
% Cada ventana de 10 s constituye una fila independiente en el Excel.
% Las tres filas de una misma transicion comparten el mismo ID_transicion.
%
% Biomarcadores por ventana:
%   RR_mean, SDNN, RMSSD, SDSD, pNN20, pNN50, CV_RR,
%   SD1, SD2, SD1_SD2,
%   DF_completo_Hz y DF_residual_Hz,
%   biomarcadores morfologicos de las ondas P y T.
%
% No se calculan LF, HF, LF/HF, LFnu, HFnu ni SampEn.
%% ============================================================

%% CONFIGURACION

% Modificar estas rutas segun la ubicacion local de las bases de datos.
bases = { ...
    struct('carpeta', ...
    'C:\Users\Emma\Documents\MATLAB\long-term-af-database-1.0.0\files', ...
    'nombre','LTAFDB'), ...
    struct('carpeta', ...
    'C:\Users\Emma\Documents\MATLAB\mit-bih-atrial-fibrillation-database-1.0.0\files', ...
    'nombre','AFDB') ...
};

carpeta_out = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_10s_RS_A_FA_3ventanas';

if ~exist(carpeta_out,'dir')
    mkdir(carpeta_out)
end

Fs_new = 500;
duracion_total_s = 180;

% Ventanas relativas al instante RS -> FA.
ventanas_rel = [ ...
    -180  -170; ...
     -90   -80; ...
     -10     0];

etiquetas_momento = [ ...
    "RS_alejado_m180_m170"; ...
    "RS_intermedio_m90_m80"; ...
    "RS_previo_m10_0"];

orden_momento = [1; 2; 3];

n_momentos = size(ventanas_rel,1);

% Controles minimos de calidad.
% En cada ventana de 10 s se exigen al menos 6 picos R, 5 intervalos RR
% validos y 6 ondas P/T validas. Las correlaciones morfologicas se
% calculan conjuntamente con todas las ondas disponibles.
min_R_180s = 60;
min_R_10s = 6;
min_RR_validos = 5;
min_ondas_10s = 6;
max_porcentaje_RR_fuera = 10;

% Se excluye porque concentra demasiadas transiciones.
registros_excluir = {'112'};

filas_resultados = {};
contador_transicion_global = 0;

%% ============================================================
% RECORRER BASES Y REGISTROS
%% ============================================================

for bb = 1:numel(bases)

    carpeta_base = bases{bb}.carpeta;
    nombre_base = bases{bb}.nombre;
    archivos_hea = dir(fullfile(carpeta_base,'*.hea'));

    fprintf('\n=============================================\n')
    fprintf('Procesando %s\n',nombre_base)
    fprintf('Registros encontrados: %d\n',numel(archivos_hea))
    fprintf('=============================================\n')

    for rr = 1:numel(archivos_hea)

        [~,nombre_registro,~] = fileparts(archivos_hea(rr).name);

        if any(strcmp(nombre_registro,registros_excluir))
            fprintf('Registro %s excluido.\n',nombre_registro)
            continue
        end

        ruta_registro = fullfile(carpeta_base,nombre_registro);
        fprintf('\n%s | %s\n',nombre_base,nombre_registro)

        try
            %% 1) LEER Y REMUESTREAR

            [sig,Fs_original] = leer_senal_wfdb_python(ruta_registro);

            if isempty(sig) || isempty(Fs_original) || ~isfinite(Fs_original)
                fprintf('  Senal no valida.\n')
                continue
            end

            x = double(sig(:,1));
            x = x(:);
            x = resample(x,Fs_new,round(Fs_original));
            Fs = Fs_new;

            %% 2) FILTRADO INICIAL
            % Se aplica el mismo criterio de filtrado y deteccion de
            % artefactos utilizado durante el preprocesamiento.
            % Las funciones auxiliares necesarias deben estar disponibles
            % en la carpeta de trabajo o en el path de MATLAB.
            % La funcion DetectarPICOSR se utiliza como funcion externa.

            [b_filt,a_filt] = butter(2,[0.5 40]/(Fs/2),'bandpass');
            x_filt = filtfilt(b_filt,a_filt,x);

            dur_artefacto = 2;

            %% 3) DETECCION DE ARTEFACTOS EN VENTANAS DE 2 s

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
            % 0 buena, 1 desconexion, 2 transicion, 3 artefacto

            %% 3.1) DETECTAR DESCONEXION

            ref_ptp_all = median(ptp_vals);
            ref_rms_all = median(rms_vals);

            low_ptp = 0.20 * ref_ptp_all;
            low_rms = 0.20 * ref_rms_all;

            idx_disc = find(ptp_vals < low_ptp | rms_vals < low_rms);
            labels(idx_disc) = 1;

            %% 3.2) REFERENCIAS DE SENAL NORMAL

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

            %% 3.3) DETECTAR TRANSICION DESPUES DE DESCONEXION

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

            %% 3.4) DETECTAR ARTEFACTO TIPO PULSOS / AMPLITUD EXTREMA

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

            %% 3.5) RELLENAR HUECOS BUENOS CORTOS ENTRE BLOQUES MALOS

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

            %% 3.6) CONVERTIR VENTANAS MALAS DE 2 s EN INTERVALOS MALOS

            intervalos_malos = [];

            for k = 1:num_vent_art

                if labels(k) ~= 0
                    intervalos_malos = [intervalos_malos; ...
                        (k-1)*dur_artefacto, k*dur_artefacto]; %#ok<AGROW>
                end
            end

            intervalos_malos = unir_intervalos(intervalos_malos);

            %% 4) SENAL PARA EL ANALISIS EN RS
            % Se utiliza la misma senal filtrada entre 0.5 y 40 Hz.
            % Dentro del calculo de DF y ondas P/T se aplica 0.5-20 Hz,
            % igual que en los analisis de RS del TFG.

            x_RS = x_filt;

            %% 5) ANOTACIONES DE RITMO

            [ann_ritmo,comments] = ...
                leer_anotaciones_wfdb_python(ruta_registro,'atr');

            if isempty(ann_ritmo) || isempty(comments)
                fprintf('  Sin anotaciones de ritmo.\n')
                continue
            end

            % Procesar las anotaciones de ritmo.
            %
            % Las anotaciones consecutivas que representan el mismo ritmo
            % se consideran parte de un unico episodio continuo.
            % Por tanto, se conserva unicamente la primera anotacion de cada
            % bloque consecutivo de RS o FA.

            ritmos_validos = {'(AFIB','(N','(NSR','(SR'};

            tiempos_ritmo = [];
            ritmos = {};

            for i = 1:min(numel(ann_ritmo),numel(comments))

                txt = strtrim(comments{i});

                if any(strcmp(txt,ritmos_validos))

                    tiempos_ritmo(end+1,1) = ...
                        double(ann_ritmo(i))/double(Fs_original); %#ok<AGROW>

                    ritmos{end+1,1} = txt; %#ok<AGROW>
                end
            end

            if numel(tiempos_ritmo) < 2
                fprintf('  Sin cambios de ritmo utilizables.\n')
                continue
            end

            [tiempos_ritmo,orden] = sort(tiempos_ritmo);
            ritmos = ritmos(orden);

            %% 5.1) UNIFICAR ANOTACIONES CONSECUTIVAS DEL MISMO RITMO
            %
            % Si aparecen varias anotaciones consecutivas que representan
            % el mismo ritmo, se consideran parte de un unico episodio.
            %
            % Ejemplo:
            %   0 s    (N
            %   400 s  (N
            %   500 s  (AFIB
            %
            % Se transforma en:
            %   0 s    (N
            %   500 s  (AFIB
            %
            % De esta forma, el inicio del episodio de RS se mantiene en
            % 0 s y no se desplaza artificialmente a 400 s.

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

            fprintf(['  Anotaciones repetidas unificadas: %d\n'], ...
                n_repetidas_eliminadas);

            fprintf(['  Episodios de ritmo tras la unificacion: %d\n'], ...
                numel(ritmos));

            if numel(tiempos_ritmo) < 2
                fprintf(['  No quedan cambios de ritmo suficientes ' ...
                    'tras unificar repeticiones.\n'])
                continue
            end

            %% 6) ANOTACIONES QRS

            ann_qrs = leer_anotaciones_muestra_wfdb_python( ...
                ruta_registro,'qrs');

            if isempty(ann_qrs)
                fprintf('  Sin anotaciones .qrs.\n')
                continue
            end

            locs_qrs_global = round(double(ann_qrs(:))*Fs/double(Fs_original));
            locs_qrs_global = unique(locs_qrs_global);
            locs_qrs_global = locs_qrs_global( ...
                locs_qrs_global >= 1 & locs_qrs_global <= numel(x_RS));

            if numel(locs_qrs_global) < 3
                continue
            end

            %% 7) BUSCAR TRANSICIONES RS -> FA

            n_validas_registro = 0;

            for ii = 1:numel(tiempos_ritmo)-1

                ritmo_actual = mapear_ritmo(ritmos{ii});
                ritmo_siguiente = mapear_ritmo(ritmos{ii+1});

                if ~(strcmp(ritmo_actual,'RS') && ...
                        strcmp(ritmo_siguiente,'FA'))
                    continue
                end

                t_ini_RS = tiempos_ritmo(ii);
                t_transicion = tiempos_ritmo(ii+1);

                % Los 180 s previos deben estar completamente dentro del
                % mismo episodio de RS.
                if (t_transicion - t_ini_RS) < duracion_total_s
                    continue
                end

                ini180 = round((t_transicion - duracion_total_s)*Fs) + 1;
                fin180 = round(t_transicion*Fs);

                if ini180 < 1 || fin180 > numel(x_RS) || fin180 <= ini180
                    continue
                end

                ecg180 = x_RS(ini180:fin180);
                ecg180 = ecg180(:);

                if any(~isfinite(ecg180))
                    continue
                end

                %% 8) COMPROBAR ARTEFACTOS SOLO EN LAS TRES VENTANAS ANALIZADAS
                %
                % No se exige que los 180 s completos esten limpios.
                % Solo se descarta la transicion si existe un artefacto en:
                %
                %   MOMENTO 1: -180 a -170 s
                %   MOMENTO 2:  -90 a  -80 s
                %   MOMENTO 3:  -10 a    0 s

                ventanas_calidad = zeros(n_momentos,2);

                for vv = 1:n_momentos
                    ventanas_calidad(vv,1) = t_transicion + ventanas_rel(vv,1);
                    ventanas_calidad(vv,2) = t_transicion + ventanas_rel(vv,2);
                end

                hay_artefacto_ventanas = false;
                ventana_con_artefacto = 0;

                for vv = 1:size(ventanas_calidad,1)

                    t_ini_comprobar = ventanas_calidad(vv,1);
                    t_fin_comprobar = ventanas_calidad(vv,2);

                    for aa = 1:size(intervalos_malos,1)

                        if intervalos_malos(aa,1) < t_fin_comprobar && ...
                                intervalos_malos(aa,2) > t_ini_comprobar

                            hay_artefacto_ventanas = true;
                            ventana_con_artefacto = vv;
                            break
                        end
                    end

                    if hay_artefacto_ventanas
                        break
                    end
                end

                if hay_artefacto_ventanas

                    fprintf(['  Transicion %.1f s omitida: la ventana %d ' ...
                        '(%g a %g s) contiene artefactos.\n'], ...
                        t_transicion, ...
                        ventana_con_artefacto, ...
                        ventanas_rel(ventana_con_artefacto,1), ...
                        ventanas_rel(ventana_con_artefacto,2));

                    continue
                end

                %% 9) R HIBRIDOS EN EL BLOQUE COMPLETO DE 180 s

                idx_qrs = locs_qrs_global >= ini180 & ...
                          locs_qrs_global <= fin180;

                qrs180 = locs_qrs_global(idx_qrs) - ini180 + 1;

                [R180,motivo_R] = DetectarPICOSR(ecg180,qrs180,Fs);

                if isempty(R180) || numel(R180) < min_R_180s
                    fprintf('  Transicion %.1f s omitida: %s\n', ...
                        t_transicion,motivo_R)
                    continue
                end

                RR180 = diff(R180)/Fs;

                if isempty(RR180)
                    continue
                end

                porcentaje_fuera = 100*mean(RR180 < 0.30 | RR180 > 1.50);

                if porcentaje_fuera > max_porcentaje_RR_fuera
                    continue
                end

                %% 10) FILTRADO ESPECIFICO DE RS SOBRE LOS 180 s

                [bRS,aRS] = butter(2,[0.5 20]/(Fs/2),'bandpass');
                ecg180_RS = filtfilt(bRS,aRS,ecg180(:));

                if isempty(ecg180_RS) || any(~isfinite(ecg180_RS))
                    fprintf(['  Transicion %.1f s omitida: ' ...
                        'ECG de 180 s filtrado no valido.\n'], ...
                        t_transicion)
                    continue
                end

                %% 11) CANCELACION QRST SOBRE LOS 180 s COMPLETOS

                residual180 = cancelar_QRST_plantilla_medianaRS( ...
                    ecg180_RS,R180,Fs);

                if isempty(residual180) || any(~isfinite(residual180))
                    fprintf(['  Transicion %.1f s omitida: ' ...
                        'residual de 180 s no valido.\n'], ...
                        t_transicion)
                    continue
                end

                residual180 = residual180(:);

                if numel(residual180) ~= numel(ecg180_RS)
                    fprintf(['  Transicion %.1f s omitida: ' ...
                        'longitud incorrecta del residual.\n'], ...
                        t_transicion)
                    continue
                end

                %% 12) TRES VENTANAS DE 10 s

                resultados_esta_transicion = cell(n_momentos,1);
                transicion_completa = true;

                for mm = 1:n_momentos

                    % Convertir tiempos relativos a indices dentro del bloque
                    % que empieza en -180 s.
                    t0_bloque = ventanas_rel(mm,1) + duracion_total_s;
                    t1_bloque = ventanas_rel(mm,2) + duracion_total_s;

                    ini10 = round(t0_bloque*Fs) + 1;
                    fin10 = round(t1_bloque*Fs);

                    if ini10 < 1 || fin10 > numel(ecg180_RS) || fin10 <= ini10
                        transicion_completa = false;
                        break
                    end

                    ecg10 = ecg180_RS(ini10:fin10);
                    ecg10 = ecg10(:);

                    residual10 = residual180(ini10:fin10);
                    residual10 = residual10(:);

                    R10 = R180(R180 >= ini10 & R180 <= fin10) - ini10 + 1;
                    R10 = unique(round(R10(:)));

                    if numel(R10) < min_R_10s
                        transicion_completa = false;
                        break
                    end

                    B = calcular_biomarcadores_RS_10s( ...
                        ecg10,residual10,R10,Fs, ...
                        min_RR_validos,min_ondas_10s);

                    if ~B.Valida
                        transicion_completa = false;
                        break
                    end

                    resultados_esta_transicion{mm} = B;
                end

                if ~transicion_completa
                    continue
                end

                contador_transicion_global = contador_transicion_global + 1;
                n_validas_registro = n_validas_registro + 1;

                id_transicion = sprintf('%s_%s_T%03d', ...
                    nombre_base,nombre_registro,n_validas_registro);

                for mm = 1:n_momentos

                    B = resultados_esta_transicion{mm};

                    filas_resultados(end+1,:) = { ... %#ok<SAGROW>
                        nombre_base,nombre_registro,id_transicion, ...
                        n_validas_registro,contador_transicion_global, ...
                        t_ini_RS,t_transicion, ...
                        etiquetas_momento(mm),orden_momento(mm), ...
                        ventanas_rel(mm,1),ventanas_rel(mm,2), ...
                        B.N_R,B.N_RR, ...
                        B.RR_mean,B.SDNN,B.RMSSD,B.SDSD, ...
                        B.pNN20,B.pNN50,B.CV_RR, ...
                        B.SD1,B.SD2,B.SD1_SD2, ...
                        B.DF_completo_Hz,B.DF_residual_Hz, ...
                        B.P_NumOndas,B.P_CorrIntraMedia,B.P_CorrIntraStd, ...
                        B.P_AmpMedia,B.P_AmpStd,B.P_StdMedia, ...
                        B.T_NumOndas,B.T_CorrIntraMedia,B.T_CorrIntraStd, ...
                        B.T_AmpMedia,B.T_AmpStd,B.T_StdMedia};
                end
            end

            fprintf('  Transiciones completas conservadas: %d\n', ...
                n_validas_registro)

        catch ME

            fprintf('  ERROR: %s\n',ME.message)

            if ~isempty(ME.stack)
                fprintf('  Linea: %d\n',ME.stack(1).line)
            end
        end
    end
end

%% ============================================================
% TABLA DE BIOMARCADORES
%% ============================================================

if isempty(filas_resultados)
    error(['No se obtuvo ninguna transicion RS->FA valida ' ...
        'con las tres ventanas de 10 s limpias.'])
end

nombres_columnas = { ...
    'Base','Registro','ID_transicion','N_transicion_registro', ...
    'N_transicion_global','Tiempo_inicio_RS_s','Tiempo_transicion_RS_FA_s', ...
    'Momento','Orden_momento','Ventana_ini_rel_s','Ventana_fin_rel_s', ...
    'N_R','N_RR', ...
    'RR_mean','SDNN','RMSSD','SDSD','pNN20','pNN50','CV_RR', ...
    'SD1','SD2','SD1_SD2','DF_completo_Hz','DF_residual_Hz', ...
    'P_NumOndas','P_CorrIntraMedia','P_CorrIntraStd', ...
    'P_AmpMedia','P_AmpStd','P_StdMedia', ...
    'T_NumOndas','T_CorrIntraMedia','T_CorrIntraStd', ...
    'T_AmpMedia','T_AmpStd','T_StdMedia'};

T = cell2table(filas_resultados,'VariableNames',nombres_columnas);

T = sortrows(T, ...
    {'Base','Registro','N_transicion_registro','Orden_momento'});

ruta_tabla = fullfile(carpeta_out, ...
    'biomarcadores_transiciones_RS_a_FA_10s_3ventanas.xlsx');

writetable(T,ruta_tabla)

fprintf('\nTabla guardada: %s\n',ruta_tabla)
fprintf('Transiciones completas: %d\n',numel(unique(T.ID_transicion)))
fprintf('Filas totales: %d\n',height(T))
fprintf('Cada transicion valida aporta tres filas de 10 s al Excel.\n')

fprintf('\nDistribucion por momento:\n')

momentos_unicos = unique(T.Momento, 'stable');

for i = 1:numel(momentos_unicos)

    momento_actual = momentos_unicos(i);

    idx = T.Momento == momento_actual;

    fprintf('  %s: %d filas\n', ...
        momento_actual, ...
        sum(idx));
end

fprintf('\nFIN.\n')

%% ============================================================
% FUNCIONES AUXILIARES PROPIAS DE ESTE SCRIPT
%% ============================================================

function B = calcular_biomarcadores_RS_10s( ...
    ecg,residual,R,Fs,min_RR,min_ondas)

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

% En una ventana de 10 s se exigen al menos 5 intervalos RR validos.

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

%% ECG DE 10 s YA FILTRADO ENTRE 0.5 Y 20 Hz

ecg_RS = ecg(:);

%% DF DEL ECG COMPLETO Y DEL RESIDUAL

[B.DF_completo_Hz,~,~,~] = ...
    frecuencia_dominante2_RS(ecg_RS,Fs);

if ~isempty(residual) && ...
        numel(residual) == numel(ecg_RS) && ...
        all(isfinite(residual))

    [B.DF_residual_Hz,~,~,~] = ...
        frecuencia_dominante2_RS(residual,Fs);
end

%% ONDA P

[ondasP,~] = extraer_ondasP_desde_R_local(ecg_RS,R,Fs);
MP = calcular_morfologia_onda_10s(ondasP,min_ondas);

B.P_NumOndas = MP.NumOndas;
B.P_CorrIntraMedia = MP.CorrIntraMedia;
B.P_CorrIntraStd = MP.CorrIntraStd;
B.P_AmpMedia = MP.AmpMedia;
B.P_AmpStd = MP.AmpStd;
B.P_StdMedia = MP.StdMedia;

%% ONDA T

[ondasT,~] = extraer_ondasT_desde_R_local(ecg_RS,R,Fs);
MT = calcular_morfologia_onda_10s(ondasT,min_ondas);

B.T_NumOndas = MT.NumOndas;
B.T_CorrIntraMedia = MT.CorrIntraMedia;
B.T_CorrIntraStd = MT.CorrIntraStd;
B.T_AmpMedia = MT.AmpMedia;
B.T_AmpStd = MT.AmpStd;
B.T_StdMedia = MT.StdMedia;

%% CONTROL FINAL
% Solo se conserva la transicion cuando todos los biomarcadores utilizados
% para el futuro UMAP son validos en las tres ventanas.

valores = [ ...
    B.RR_mean B.SDNN B.RMSSD B.SDSD B.pNN20 B.pNN50 ...
    B.CV_RR B.SD1 B.SD2 B.SD1_SD2 ...
    B.DF_completo_Hz B.DF_residual_Hz ...
    B.P_CorrIntraMedia B.P_CorrIntraStd ...
    B.P_AmpMedia B.P_AmpStd B.P_StdMedia ...
    B.T_CorrIntraMedia B.T_CorrIntraStd ...
    B.T_AmpMedia B.T_AmpStd B.T_StdMedia];

B.Valida = all(isfinite(valores));
end

function M = calcular_morfologia_onda_10s(ondas,min_ondas)

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

function ritmo = mapear_ritmo(txt)
    if strcmp(txt,'(AFIB')
        ritmo = 'FA';
    elseif strcmp(txt,'(N') || strcmp(txt,'(NSR') || strcmp(txt,'(SR')
        ritmo = 'RS';
    else
        ritmo = 'OTHER';
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
