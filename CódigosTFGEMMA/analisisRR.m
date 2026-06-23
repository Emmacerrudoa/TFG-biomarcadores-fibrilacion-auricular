clear
clc
close all

%% ============================================================
% CÁLCULO DE BIOMARCADORES TEMPORALES A PARTIR DE INTERVALOS RR
%
% Este script calcula biomarcadores de variabilidad del ritmo cardíaco
% a partir de las ventanas ECG generadas previamente.
%
% Para cada grupo se recorren los archivos .mat, se agrupan las ventanas
% por registro/paciente y se calculan biomarcadores globales y biomarcadores
% dependientes de continuidad temporal.
%
% Biomarcadores calculados:
%   - RR_mean, SDNN y CV_RR.
%   - RMSSD, SDSD, pNN50 y pNN20.
%   - SD1, SD2 y SD1/SD2 mediante análisis de Poincaré.
%   - SampEn.
%   - LF, HF, LF/HF, LFnu y HFnu.
%
% Salidas:
%   - resultados_RR_4grupos.mat
%   - biomarcadores_RR_FINAL_con_poincare_sampen_frecuencia.xlsx
%
% Requisitos:
%   - Archivos .mat generados en la fase de segmentación.
%   - MATLAB con Signal Processing Toolbox.
%% ============================================================

%% CONFIGURACIÓN GENERAL

Fs_default = 500;

% Modificar estas rutas según la ubicación local del dataset generado.
carpetas = { ...
    struct('ruta','C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\SANO','grupo','SANO'), ...
    struct('ruta','C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\FA_PAROXISTICA_RS','grupo','FA_PAROXISTICA_RS'), ...
    struct('ruta','C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\FA_PAROXISTICA_FA','grupo','FA_PAROXISTICA_FA'), ...
    struct('ruta','C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO\FA_PERSISTENTE','grupo','FA_PERSISTENTE') ...
};

carpeta_out = 'C:\Users\Emma\Documents\MATLAB\analisisRR_hibrido';
if ~exist(carpeta_out,'dir')
    mkdir(carpeta_out);
end

edges_rr = 0.30:0.05:1.50;
tol_cont = 1/Fs_default;

%% ESTRUCTURAS DE RESULTADOS

resultados_por_grupo = struct([]);

%% 1. PROCESAR CADA GRUPO

for g = 1:numel(carpetas)

    ruta_grupo = carpetas{g}.ruta;
    nombre_grupo = carpetas{g}.grupo;

    fprintf('\n=============================\n');
    fprintf('Procesando grupo: %s\n', nombre_grupo);
    fprintf('=============================\n');

    archivos = dir(fullfile(ruta_grupo, '*.mat'));

    if isempty(archivos)
        fprintf('No se encontraron archivos .mat en %s\n', ruta_grupo);
        resultados_por_grupo(g).grupo = nombre_grupo;
        resultados_por_grupo(g).resultados = struct([]);
        continue
    end

    %% Leer nombres de los registros

    registros = cell(numel(archivos),1);
    for k = 1:numel(archivos)
        D = load(fullfile(ruta_grupo, archivos(k).name), 'nombre_registro');

        if isfield(D,'nombre_registro') && ~isempty(D.nombre_registro)
            registros{k} = D.nombre_registro;
        else
            registros{k} = '';
        end
    end

    registros_unicos = unique(registros);
    registros_unicos(cellfun(@isempty, registros_unicos)) = [];

    resultados_grupo = struct([]);

    %% 2. RECORRER PACIENTES / REGISTROS

    for r = 1:numel(registros_unicos)

        reg = registros_unicos{r};
        fprintf('Procesando %s - %s\n', nombre_grupo, reg);

        idx_reg = strcmp(registros, reg);
        archivos_reg = archivos(idx_reg);

        %% 2.1 CARGAR INFORMACIÓN DE LAS VENTANAS

        info_ventanas = struct([]);
        cont_info = 0;

        for k = 1:numel(archivos_reg)

            D = load(fullfile(ruta_grupo, archivos_reg(k).name), ...
                'locs_R','t_ini','t_fin','nombre_registro','Fs','nombre_base', ...
                'tipo_registro','ritmo_ventana','ID_global');

            if ~isfield(D,'locs_R') || ~isfield(D,'t_ini') || ~isfield(D,'t_fin')
                continue
            end

            if isempty(D.locs_R) || ~isfinite(D.t_ini) || ~isfinite(D.t_fin)
                continue
            end

            cont_info = cont_info + 1;

            info_ventanas(cont_info).archivo = archivos_reg(k).name;
            info_ventanas(cont_info).ruta = fullfile(ruta_grupo, archivos_reg(k).name);
            info_ventanas(cont_info).locs_R = D.locs_R(:);
            info_ventanas(cont_info).t_ini = D.t_ini;
            info_ventanas(cont_info).t_fin = D.t_fin;

            if isfield(D,'Fs') && ~isempty(D.Fs) && isfinite(D.Fs)
                info_ventanas(cont_info).Fs = D.Fs;
            else
                info_ventanas(cont_info).Fs = Fs_default;
            end

            if isfield(D,'nombre_base')
                info_ventanas(cont_info).nombre_base = D.nombre_base;
            else
                info_ventanas(cont_info).nombre_base = '';
            end

            if isfield(D,'tipo_registro')
                info_ventanas(cont_info).tipo_registro = D.tipo_registro;
            else
                info_ventanas(cont_info).tipo_registro = '';
            end

            if isfield(D,'ritmo_ventana')
                info_ventanas(cont_info).ritmo_ventana = D.ritmo_ventana;
            else
                info_ventanas(cont_info).ritmo_ventana = '';
            end

            if isfield(D,'ID_global')
                info_ventanas(cont_info).ID_global = D.ID_global;
            else
                info_ventanas(cont_info).ID_global = NaN;
            end
        end

        if isempty(info_ventanas)
            fprintf('No hay ventanas validas en %s\n', reg);
            continue
        end

        %% 2.2 ORDENAR POR TIEMPO DE INICIO

        [~, ord] = sort([info_ventanas.t_ini]);
        info_ventanas = info_ventanas(ord);

        %% ============================================================
        % RAMA A: RR GLOBALES POR VENTANA
        % RR_mean, SDNN, CV_RR e histograma RR.
        % ============================================================

        RR_total_global = [];
        n_ventanas_validas_global = 0;

        %% ============================================================
        % RAMA B: BLOQUES CONTINUOS EN TIEMPO ABSOLUTO
        % RMSSD, SDSD, pNN50, pNN20, SD1, SD2, SampEn, LF, HF y LF/HF.
        % ============================================================

        RR_bloques = {};
        R_abs_bloque = [];
        prev_t_fin = NaN;
        primer_bloque = true;

        for k = 1:numel(info_ventanas)

            locs_R = info_ventanas(k).locs_R(:);
            t_ini = info_ventanas(k).t_ini;
            t_fin = info_ventanas(k).t_fin;
            Fs = info_ventanas(k).Fs;

            if isempty(locs_R) || numel(locs_R) < 3 || ~isfinite(Fs) || Fs <= 0
                continue
            end

            %% RAMA A: RR por ventana

            RR_vent = diff(locs_R) / Fs;
            RR_vent = RR_vent(isfinite(RR_vent));
            RR_vent = RR_vent(RR_vent >= 0.30 & RR_vent <= 1.50);

            if ~isempty(RR_vent)
                RR_total_global = [RR_total_global; RR_vent(:)]; %#ok<AGROW>
                n_ventanas_validas_global = n_ventanas_validas_global + 1;
            end

            %% RAMA B: tiempos absolutos por bloque

            R_abs = t_ini + (locs_R(:)-1)/Fs;
            R_abs = R_abs(isfinite(R_abs));

            if isempty(R_abs)
                continue
            end

            if primer_bloque
                R_abs_bloque = R_abs;
                primer_bloque = false;
            else
                es_continua = abs(t_ini - prev_t_fin) <= tol_cont;

                if es_continua
                    R_abs_bloque = [R_abs_bloque; R_abs]; %#ok<AGROW>
                else
                    RR_bloque = diff(sort(R_abs_bloque));
                    RR_bloque = RR_bloque(isfinite(RR_bloque));
                    RR_bloque = RR_bloque(RR_bloque >= 0.30 & RR_bloque <= 1.50);

                    if ~isempty(RR_bloque)
                        RR_bloques{end+1,1} = RR_bloque; %#ok<SAGROW>
                    end

                    R_abs_bloque = R_abs;
                end
            end

            prev_t_fin = t_fin;
        end

        %% Cerrar último bloque continuo

        if ~primer_bloque
            RR_bloque = diff(sort(R_abs_bloque));
            RR_bloque = RR_bloque(isfinite(RR_bloque));
            RR_bloque = RR_bloque(RR_bloque >= 0.30 & RR_bloque <= 1.50);

            if ~isempty(RR_bloque)
                RR_bloques{end+1,1} = RR_bloque;
            end
        end

        %% 2.3 COMPROBAR SI HAY SUFICIENTES INTERVALOS RR

        if numel(RR_total_global) < 50
            fprintf('Muy pocos RR validos en %s\n', reg);
            continue
        end

        %% 3. BIOMARCADORES RR GLOBALES POR PACIENTE

        RR_mean = mean(RR_total_global, 'omitnan');
        SDNN = std(RR_total_global, 0, 'omitnan');

        if RR_mean == 0 || isnan(RR_mean)
            CV_RR = NaN;
        else
            CV_RR = SDNN / RR_mean;
        end

        %% 4. BIOMARCADORES CALCULADOS EN BLOQUES CONTINUOS

        dRR = [];
        RRn = [];
        RRn1 = [];
        n_segmentos_validos = 0;
        segmentos_RR = {};

        SampEn_vals = [];
        LF_vals = [];
        HF_vals = [];
        LF_HF_vals = [];
        LFnu_vals = [];
        HFnu_vals = [];

        for b = 1:numel(RR_bloques)

            RRb = RR_bloques{b};
            RRb = RRb(:);
            RRb = RRb(isfinite(RRb));
            RRb = RRb(RRb >= 0.30 & RRb <= 1.50);

            if isempty(RRb)
                continue
            end

            segmentos_RR{end+1,1} = RRb(:); %#ok<SAGROW>
            n_segmentos_validos = n_segmentos_validos + 1;

            %% Diferencias consecutivas y análisis de Poincaré

            if numel(RRb) >= 2
                dRR = [dRR; diff(RRb(:))]; %#ok<AGROW>
                RRn = [RRn; RRb(1:end-1)]; %#ok<AGROW>
                RRn1 = [RRn1; RRb(2:end)]; %#ok<AGROW>
            end

            %% Sample Entropy por bloque continuo

            RRb_sampen = RRb;

            if numel(RRb_sampen) > 1000
                ini_sampen = floor((numel(RRb_sampen) - 1000)/2) + 1;
                RRb_sampen = RRb_sampen(ini_sampen:ini_sampen+999);
            end
            
            se_b = calcular_sampen_rr(RRb_sampen, 2, 0.20);

            if isfinite(se_b)
                SampEn_vals(end+1,1) = se_b; %#ok<AGROW>
            end

            %% Análisis en frecuencia por bloque continuo

            [LF_b, HF_b, LF_HF_b, LFnu_b, HFnu_b] = calcular_lfhf_rr(RRb);
            
            if isfinite(LF_b)
                LF_vals(end+1,1) = LF_b; %#ok<AGROW>
            end
            
            if isfinite(HF_b)
                HF_vals(end+1,1) = HF_b; %#ok<AGROW>
            end
            
            if isfinite(LF_HF_b)
                LF_HF_vals(end+1,1) = LF_HF_b; %#ok<AGROW>
            end
            
            if isfinite(LFnu_b)
                LFnu_vals(end+1,1) = LFnu_b; %#ok<AGROW>
            end
            
            if isfinite(HFnu_b)
                HFnu_vals(end+1,1) = HFnu_b; %#ok<AGROW>
            end
        end

        %% 4.1 RMSSD, SDSD, pNN50, pNN20

        if isempty(dRR)
            RMSSD = NaN;
            SDSD = NaN;
            pNN50 = NaN;
            pNN20 = NaN;
        else
            RMSSD = sqrt(mean(dRR.^2, 'omitnan'));
            SDSD  = std(dRR, 0, 'omitnan');
            pNN50 = 100 * sum(abs(dRR) > 0.05) / numel(dRR);
            pNN20 = 100 * sum(abs(dRR) > 0.02) / numel(dRR);
        end

        %% 4.2 Sample Entropy

        if isempty(SampEn_vals)
            SampEn = NaN;
        else
            SampEn = mean(SampEn_vals, 'omitnan');
        end

        %% 4.3 Análisis en frecuencia

        if isempty(LF_vals)
            LF = NaN;
        else
            LF = mean(LF_vals, 'omitnan');
        end

        if isempty(HF_vals)
            HF = NaN;
        else
            HF = mean(HF_vals, 'omitnan');
        end

        if isempty(LF_HF_vals)
            LF_HF = NaN;
        else
            LF_HF = mean(LF_HF_vals, 'omitnan');
        end

        if isempty(LFnu_vals)
            LFnu = NaN;
        else
            LFnu = mean(LFnu_vals, 'omitnan');
        end
        
        if isempty(HFnu_vals)
            HFnu = NaN;
        else
            HFnu = mean(HFnu_vals, 'omitnan');
        end
       
        %% 5. ANÁLISIS DE POINCARÉ CON PARES INTRA-BLOQUE

        if isempty(RRn)
            SD1 = NaN;
            SD2 = NaN;
            SD1_SD2_ratio = NaN;
        else
            xp = (RRn + RRn1) / sqrt(2);
            yp = (RRn1 - RRn) / sqrt(2);

            SD1 = std(yp, 0, 'omitnan');
            SD2 = std(xp, 0, 'omitnan');

            if SD2 == 0 || isnan(SD2)
                SD1_SD2_ratio = NaN;
            else
                SD1_SD2_ratio = SD1 / SD2;
            end
        end

        %% 6. HISTOGRAMA RR

        hist_rr = histcounts(RR_total_global, edges_rr, 'Normalization', 'probability');

        %% 7. GUARDAR RESULTADOS DEL PACIENTE

        idx_res = numel(resultados_grupo) + 1;

        resultados_grupo(idx_res).registro = reg;
        resultados_grupo(idx_res).grupo = nombre_grupo;

        resultados_grupo(idx_res).RR = RR_total_global;
        resultados_grupo(idx_res).RR_mean = RR_mean;
        resultados_grupo(idx_res).SDNN = SDNN;
        resultados_grupo(idx_res).CV_RR = CV_RR;
        resultados_grupo(idx_res).hist_rr = hist_rr;
        resultados_grupo(idx_res).n_RR = numel(RR_total_global);
        resultados_grupo(idx_res).n_ventanas_validas = n_ventanas_validas_global;

        resultados_grupo(idx_res).RMSSD = RMSSD;
        resultados_grupo(idx_res).SDSD = SDSD;
        resultados_grupo(idx_res).pNN50 = pNN50;
        resultados_grupo(idx_res).pNN20 = pNN20;
        resultados_grupo(idx_res).n_dRR = numel(dRR);

        resultados_grupo(idx_res).SD1 = SD1;
        resultados_grupo(idx_res).SD2 = SD2;
        resultados_grupo(idx_res).SD1_SD2_ratio = SD1_SD2_ratio;

        resultados_grupo(idx_res).SampEn = SampEn;
        resultados_grupo(idx_res).LF = LF;
        resultados_grupo(idx_res).HF = HF;
        resultados_grupo(idx_res).LF_HF = LF_HF;
        resultados_grupo(idx_res).LFnu = LFnu;
        resultados_grupo(idx_res).HFnu = HFnu;

        resultados_grupo(idx_res).n_segmentos = n_segmentos_validos;
        resultados_grupo(idx_res).segmentos_RR = segmentos_RR;
    end

    %% GUARDAR RESULTADOS DEL GRUPO

    resultados_por_grupo(g).grupo = nombre_grupo;
    resultados_por_grupo(g).resultados = resultados_grupo;
end

%% 8. GUARDAR .MAT GLOBAL

save(fullfile(carpeta_out, 'resultados_RR_4grupos.mat'), ...
    'resultados_por_grupo', 'edges_rr', '-v7.3');

%% 9. TABLA FINAL

registro_tabla = {};
grupo_tabla = {};
n_rr_tabla = [];
n_ventanas_validas_tabla = [];
n_segmentos_tabla = [];
n_drr_tabla = [];

rr_mean_tabla = [];
sdnn_tabla = [];
rmssd_tabla = [];
sdsd_tabla = [];
pnn50_tabla = [];
pnn20_tabla = [];
cv_rr_tabla = [];
sd1_tabla = [];
sd2_tabla = [];
sd1sd2_ratio_tabla = [];

sampen_tabla = [];
lf_tabla = [];
hf_tabla = [];
lfhf_tabla = [];
lfnu_tabla = [];
hfnu_tabla = [];

for g = 1:numel(resultados_por_grupo)
    Rg = resultados_por_grupo(g).resultados;
    nombre_grupo = resultados_por_grupo(g).grupo;

    for i = 1:numel(Rg)
        registro_tabla{end+1,1} = Rg(i).registro; %#ok<SAGROW>
        grupo_tabla{end+1,1} = nombre_grupo; %#ok<SAGROW>
        n_rr_tabla(end+1,1) = Rg(i).n_RR; %#ok<SAGROW>
        n_ventanas_validas_tabla(end+1,1) = Rg(i).n_ventanas_validas; %#ok<SAGROW>
        n_segmentos_tabla(end+1,1) = Rg(i).n_segmentos; %#ok<SAGROW>
        n_drr_tabla(end+1,1) = Rg(i).n_dRR; %#ok<SAGROW>

        rr_mean_tabla(end+1,1) = Rg(i).RR_mean; %#ok<SAGROW>
        sdnn_tabla(end+1,1) = Rg(i).SDNN; %#ok<SAGROW>
        rmssd_tabla(end+1,1) = Rg(i).RMSSD; %#ok<SAGROW>
        sdsd_tabla(end+1,1) = Rg(i).SDSD; %#ok<SAGROW>
        pnn50_tabla(end+1,1) = Rg(i).pNN50; %#ok<SAGROW>
        pnn20_tabla(end+1,1) = Rg(i).pNN20; %#ok<SAGROW>
        cv_rr_tabla(end+1,1) = Rg(i).CV_RR; %#ok<SAGROW>
        sd1_tabla(end+1,1) = Rg(i).SD1; %#ok<SAGROW>
        sd2_tabla(end+1,1) = Rg(i).SD2; %#ok<SAGROW>
        sd1sd2_ratio_tabla(end+1,1) = Rg(i).SD1_SD2_ratio; %#ok<SAGROW>

        sampen_tabla(end+1,1) = Rg(i).SampEn; %#ok<SAGROW>
        lf_tabla(end+1,1) = Rg(i).LF; %#ok<SAGROW>
        hf_tabla(end+1,1) = Rg(i).HF; %#ok<SAGROW>
        lfhf_tabla(end+1,1) = Rg(i).LF_HF; %#ok<SAGROW>
        lfnu_tabla(end+1,1) = Rg(i).LFnu; %#ok<SAGROW>
        hfnu_tabla(end+1,1) = Rg(i).HFnu; %#ok<SAGROW>
    end
end

if ~isempty(registro_tabla)
    T_bio = table(registro_tabla, grupo_tabla, n_rr_tabla, ...
        n_ventanas_validas_tabla, n_segmentos_tabla, n_drr_tabla, ...
        rr_mean_tabla, sdnn_tabla, rmssd_tabla, sdsd_tabla, ...
        pnn50_tabla, pnn20_tabla, cv_rr_tabla, ...
        sd1_tabla, sd2_tabla, sd1sd2_ratio_tabla, ...
        sampen_tabla, lf_tabla, hf_tabla, lfhf_tabla, lfnu_tabla, hfnu_tabla, ...
        'VariableNames', {'Paciente','Grupo','NumeroRR','NumeroVentanasValidas', ...
        'NumeroSegmentos','Numero_dRR','RR_mean','SDNN','RMSSD','SDSD', ...
        'pNN50','pNN20','CV_RR','SD1','SD2','SD1_SD2_ratio', ...
        'SampEn','LF','HF','LF_HF','LFnu','HFnu'});

    writetable(T_bio, fullfile(carpeta_out, 'biomarcadores_RR_FINAL_con_poincare_sampen_frecuencia.xlsx'));
end

fprintf('\nProceso terminado.\n');

function [LF, HF, LF_HF, LFnu, HFnu, f, Pxx] = calcular_lfhf_rr(RR)

LF = NaN;
HF = NaN;
LF_HF = NaN;
LFnu = NaN;
HFnu = NaN;
f = [];
Pxx = [];

RR = RR(:);
RR = RR(isfinite(RR));
RR = RR(RR >= 0.30 & RR <= 1.50);

if numel(RR) < 120
    return
end

t_rr = cumsum(RR);
t_rr = t_rr - t_rr(1);

if t_rr(end) < 300
    return
end

fs_interp = 4;
t_uniforme = (t_rr(1):1/fs_interp:t_rr(end))';

rr_interp = interp1(t_rr, RR, t_uniforme, 'linear');

idx_ok = isfinite(rr_interp);
if sum(idx_ok) < 256
    return
end

rr_interp = rr_interp(idx_ok);

rr_interp = detrend(rr_interp, 'linear');
rr_interp = rr_interp - mean(rr_interp, 'omitnan');

if std(rr_interp, 0, 'omitnan') <= 0
    return
end

win = 512;

noverlap = round(0.5 * win);

[Pxx, f] = pwelch(rr_interp, hamming(win), noverlap, [], fs_interp);

idx_LF = (f >= 0.04) & (f < 0.15);
idx_HF = (f >= 0.15) & (f <= 0.40);

if ~any(idx_LF) || ~any(idx_HF)
    return
end

LF = trapz(f(idx_LF), Pxx(idx_LF));
HF = trapz(f(idx_HF), Pxx(idx_HF));

if isfinite(HF) && HF > 0
    LF_HF = LF / HF;
end

TP = LF + HF;

if isfinite(TP) && TP > 0
    LFnu = 100 * LF / TP;
    HFnu = 100 * HF / TP;
end

end

function sampen = calcular_sampen_rr(RR, m, r_factor)

% Calcula la Sample Entropy (SampEn) de una serie RR.
%
% RR       -> serie de intervalos RR en segundos
% m        -> dimensión del patrón
% r_factor -> tolerancia relativa a la desviación estándar de RR

sampen = NaN;

RR = RR(:);
RR = RR(isfinite(RR));

if numel(RR) < 20
    return
end

sdRR = std(RR, 0, 'omitnan');
if ~isfinite(sdRR) || sdRR <= 0
    return
end

r = r_factor * sdRR;
N = length(RR);

A = 0;
B = 0;

for i = 1:(N - m)
    xi_m = RR(i:i+m-1);
    xi_m1 = RR(i:i+m);

    for j = (i+1):(N - m)
        xj_m = RR(j:j+m-1);

        if max(abs(xi_m - xj_m)) <= r
            B = B + 1;

            xj_m1 = RR(j:j+m);
            if max(abs(xi_m1 - xj_m1)) <= r
                A = A + 1;
            end
        end
    end
end

if B == 0 || A == 0
    return
end

sampen = -log(A / B);

end