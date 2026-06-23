clear
clc
close all

%% ============================================================
% FRECUENCIA DOMINANTE EN ECG COMPLETO Y RESIDUAL
%
% Este script calcula la frecuencia dominante (DF) en ventanas de 120 s
% del dataset final. El análisis se realiza tanto sobre el ECG completo
% como sobre la señal residual obtenida tras la cancelación QRS-T mediante
% una plantilla mediana.
%
% La DF se calcula inicialmente en cada ventana y después se agrupa por
% paciente mediante la combinación Base + Registro. El resumen final por
% grupo se obtiene a partir de la media de cada paciente, de forma que cada
% paciente aporta un único valor al análisis grupal.
%
% Grupos analizados:
%   - SANO
%   - FA_PAROXISTICA_RS
%   - FA_PAROXISTICA_FA
%   - FA_PERSISTENTE
%
% Salidas:
%   - resultados_DF_por_ventana_120s.xlsx
%   - resultados_DF_por_paciente.xlsx
%   - resumen_DF_por_grupo_pacientes.xlsx
%   - boxplots de DF por paciente
%   - histogramas de DF por grupo
%
% Funciones locales incluidas:
%   - limpiar_locs_local
%   - frecuencia_dominante2_RS120s
%   - frecuencia_dominante2_FA120s
%   - cancelar_QRST_plantilla_medianaRS
%   - cancelar_QRST_plantilla_medianaFA
%
% Requisitos:
%   - Dataset final en formato .mat.
%   - Signal Processing Toolbox.
%% ============================================================

%% CONFIGURACIÓN GENERAL

carpeta_dataset = ...
    'C:\Users\Emma\Documents\MATLAB\dataset_final_HIBRIDO';

carpeta_out = ...
    'C:\Users\Emma\Documents\MATLAB\DF_120s_completo_y_residual';

if ~exist(carpeta_out, 'dir')
    mkdir(carpeta_out);
end

grupos = { ...
    'FA_PAROXISTICA_FA', ...
    'FA_PAROXISTICA_RS', ...
    'SANO', ...
    'FA_PERSISTENTE'};

min_R = 15;

resultados = {};

%% ============================================================
% RECORRER LOS GRUPOS DEL DATASET
%% ============================================================

for g = 1:numel(grupos)

    grupo = grupos{g};

    carpeta_grupo = fullfile( ...
        carpeta_dataset, ...
        grupo);

    if ~exist(carpeta_grupo, 'dir')

        warning('No existe la carpeta: %s', carpeta_grupo);
        continue
    end

    archivos = dir(fullfile( ...
        carpeta_grupo, ...
        '*.mat'));

    fprintf('\n============================================\n');
    fprintf('Procesando grupo: %s\n', grupo);
    fprintf('Ventanas encontradas: %d\n', numel(archivos));
    fprintf('============================================\n');

    %% Definir ritmo según el grupo

    if strcmp(grupo, 'SANO') || ...
            strcmp(grupo, 'FA_PAROXISTICA_RS')

        ritmo = "RS";
        es_RS = true;
        es_FA = false;

    elseif strcmp(grupo, 'FA_PAROXISTICA_FA') || ...
            strcmp(grupo, 'FA_PERSISTENTE')

        ritmo = "FA";
        es_RS = false;
        es_FA = true;

    else

        fprintf('Grupo no reconocido: %s\n', grupo);
        continue
    end

    %% Recorrer ventanas del grupo

    for k = 1:numel(archivos)

        if mod(k,50) == 0 || k == 1

            fprintf( ...
                '  %s: ventana %d/%d | archivo: %s\n', ...
                grupo, ...
                k, ...
                numel(archivos), ...
                archivos(k).name);
        end

        archivo_mat = fullfile( ...
            carpeta_grupo, ...
            archivos(k).name);

        try

            S = load(archivo_mat);

            %% ------------------------------------------------------------
            % Leer señal
            %% ------------------------------------------------------------

            if ~isfield(S, 'ventana')

                fprintf( ...
                    '  Omitida %s: no contiene variable ventana.\n', ...
                    archivos(k).name);

                continue
            end

            x = double(S.ventana(:));

            if isempty(x) || any(~isfinite(x))

                fprintf( ...
                    '  Omitida %s: señal vacía o no finita.\n', ...
                    archivos(k).name);

                continue
            end

            %% ------------------------------------------------------------
            % Leer frecuencia de muestreo
            %% ------------------------------------------------------------

            if isfield(S, 'Fs')
                Fs = double(S.Fs);
            else
                Fs = 500;
            end

            if isempty(Fs) || ~isfinite(Fs) || Fs <= 0

                fprintf( ...
                    '  Omitida %s: Fs no válido.\n', ...
                    archivos(k).name);

                continue
            end

            %% ------------------------------------------------------------
            % Leer picos R
            %% ------------------------------------------------------------

            if ~isfield(S, 'locs_R')

                fprintf( ...
                    '  Omitida %s: no contiene locs_R.\n', ...
                    archivos(k).name);

                continue
            end

            locs_R = limpiar_locs_local( ...
                S.locs_R, ...
                length(x));

            if numel(locs_R) < min_R

                fprintf( ...
                    '  Omitida %s: pocos picos R.\n', ...
                    archivos(k).name);

                continue
            end

            %% ------------------------------------------------------------
            % Metadatos
            %% ------------------------------------------------------------

            nombre_registro = "";
            nombre_base = "";
            tipo_registro = "";
            ritmo_ventana = "";

            ID_global = NaN;
            t_ini = NaN;
            t_fin = NaN;
            s = NaN;
            w = NaN;

            if isfield(S, 'nombre_registro')
                nombre_registro = string(S.nombre_registro);
            end

            if isfield(S, 'nombre_base')
                nombre_base = string(S.nombre_base);
            end

            if isfield(S, 'tipo_registro')
                tipo_registro = string(S.tipo_registro);
            end

            if isfield(S, 'ritmo_ventana')
                ritmo_ventana = string(S.ritmo_ventana);
            end

            if isfield(S, 'ID_global')
                ID_global = double(S.ID_global);
            end

            if isfield(S, 't_ini')
                t_ini = double(S.t_ini);
            end

            if isfield(S, 't_fin')
                t_fin = double(S.t_fin);
            end

            if isfield(S, 's')
                s = double(S.s);
            end

            if isfield(S, 'w')
                w = double(S.w);
            end

            %% ------------------------------------------------------------
            % Comprobar que Base y Registro estén disponibles
            %% ------------------------------------------------------------

            if strlength(nombre_base) == 0 || ...
                    strlength(nombre_registro) == 0

                fprintf( ...
                    '  Omitida %s: falta Base o Registro.\n', ...
                    archivos(k).name);

                continue
            end

            %% ------------------------------------------------------------
            % ECG completo para análisis
            %
            % RS:
            %   se filtra entre 0.5 y 20 Hz.
            %
            % FA:
            %   se utiliza la señal guardada.
            %% ------------------------------------------------------------

            if es_RS

                [b_RS, a_RS] = butter( ...
                    2, ...
                    [0.5 20] / (Fs/2), ...
                    'bandpass');

                x_completo = filtfilt( ...
                    b_RS, ...
                    a_RS, ...
                    x);

            elseif es_FA

                x_completo = x;
            end

            x_completo = double(x_completo(:));

            %% ------------------------------------------------------------
            % 1) DF DEL ECG COMPLETO
            %% ------------------------------------------------------------

            DF_completo = NaN;
            Power_completo = NaN;

            try

                if es_RS

                    [DF_completo, ~, ~, Power_completo] = ...
                        frecuencia_dominante2_RS120s( ...
                        x_completo, ...
                        Fs);

                elseif es_FA

                    [DF_completo, ~, ~, Power_completo] = ...
                        frecuencia_dominante2_FA120s( ...
                        x_completo, ...
                        Fs);
                end

            catch ME_DF_comp

                fprintf( ...
                    '  ERROR DF completo: %s | %s\n', ...
                    archivos(k).name, ...
                    ME_DF_comp.message);

                DF_completo = NaN;
                Power_completo = NaN;
            end

            %% ------------------------------------------------------------
            % 2) CANCELACIÓN QRS-T Y DF RESIDUAL
            %% ------------------------------------------------------------

            DF_residual = NaN;
            Power_residual = NaN;
            residual = [];

            try

                if es_RS

                    residual = ...
                        cancelar_QRST_plantilla_medianaRS( ...
                        x_completo, ...
                        locs_R, ...
                        Fs);

                elseif es_FA

                    residual = ...
                        cancelar_QRST_plantilla_medianaFA( ...
                        x_completo, ...
                        locs_R, ...
                        Fs);
                end

                if ~isempty(residual)

                    residual = double(residual(:));

                    if all(isfinite(residual))

                        if es_RS

                            [DF_residual, ~, ~, Power_residual] = ...
                                frecuencia_dominante2_RS120s( ...
                                residual, ...
                                Fs);

                        elseif es_FA

                            [DF_residual, ~, ~, Power_residual] = ...
                                frecuencia_dominante2_FA120s( ...
                                residual, ...
                                Fs);
                        end
                    end
                end

            catch ME_cancel

                fprintf( ...
                    '  ERROR cancelando QRST: %s | %s\n', ...
                    archivos(k).name, ...
                    ME_cancel.message);

                DF_residual = NaN;
                Power_residual = NaN;
            end

            %% ------------------------------------------------------------
            % Guardar resultado de la ventana
            %% ------------------------------------------------------------

            resultados(end+1,:) = { ...
                string(grupo), ...
                ritmo, ...
                string(archivos(k).name), ...
                nombre_base, ...
                nombre_registro, ...
                tipo_registro, ...
                ritmo_ventana, ...
                ID_global, ...
                s, ...
                w, ...
                t_ini, ...
                t_fin, ...
                Fs, ...
                length(x) / Fs, ...
                numel(locs_R), ...
                DF_completo, ...
                Power_completo, ...
                DF_residual, ...
                Power_residual}; %#ok<SAGROW>

        catch ME

            fprintf( ...
                '  ERROR en %s: %s\n', ...
                archivos(k).name, ...
                ME.message);

            if ~isempty(ME.stack)

                fprintf( ...
                    '  Línea aproximada: %d\n', ...
                    ME.stack(1).line);
            end
        end
    end
end

%% ============================================================
% TABLA DE RESULTADOS POR VENTANA
%% ============================================================

if isempty(resultados)

    warning('No se generaron resultados.');
    return
end

T = cell2table( ...
    resultados, ...
    'VariableNames', { ...
    'Grupo', ...
    'Ritmo', ...
    'Archivo', ...
    'Base', ...
    'Registro', ...
    'TipoRegistro', ...
    'RitmoVentana', ...
    'ID_global', ...
    'Segmento_s', ...
    'Ventana_w', ...
    't_ini_s', ...
    't_fin_s', ...
    'Fs', ...
    'Duracion_s', ...
    'N_R', ...
    'DF_completo_Hz', ...
    'Power_completo', ...
    'DF_residual_Hz', ...
    'Power_residual'});

T.Grupo = string(T.Grupo);
T.Ritmo = string(T.Ritmo);
T.Archivo = string(T.Archivo);
T.Base = string(T.Base);
T.Registro = string(T.Registro);
T.TipoRegistro = string(T.TipoRegistro);
T.RitmoVentana = string(T.RitmoVentana);

%% Crear identificador de paciente mediante Base + Registro

T.Paciente = T.Base + "_" + T.Registro;

writetable( ...
    T, ...
    fullfile( ...
    carpeta_out, ...
    'resultados_DF_por_ventana_120s.xlsx'));

%% ============================================================
% AGRUPAR RESULTADOS POR PACIENTE
%
% Ejemplos:
%
%   BASE_1 + 00 = BASE_1_00
%   BASE_2 + 00 = BASE_2_00
%
% Aunque el número de registro sea igual, son pacientes distintos
% porque pertenecen a bases diferentes.
%% ============================================================

claves = strcat( ...
    T.Grupo, ...
    "__", ...
    T.Paciente);

[grupo_paciente, claves_unicas] = findgroups(claves);

num_pacientes = numel(claves_unicas);

Grupo = strings(num_pacientes, 1);
Ritmo = strings(num_pacientes, 1);
Paciente = strings(num_pacientes, 1);
Base = strings(num_pacientes, 1);
Registro = strings(num_pacientes, 1);

N_ventanas_total = zeros(num_pacientes, 1);

N_DF_completo = zeros(num_pacientes, 1);
DF_completo_media_Hz = NaN(num_pacientes, 1);
DF_completo_std_Hz = NaN(num_pacientes, 1);

N_DF_residual = zeros(num_pacientes, 1);
DF_residual_media_Hz = NaN(num_pacientes, 1);
DF_residual_std_Hz = NaN(num_pacientes, 1);

for i = 1:num_pacientes

    idx = grupo_paciente == i;

    primera_fila = find(idx, 1);

    Grupo(i) = T.Grupo(primera_fila);
    Ritmo(i) = T.Ritmo(primera_fila);
    Paciente(i) = T.Paciente(primera_fila);
    Base(i) = T.Base(primera_fila);
    Registro(i) = T.Registro(primera_fila);

    N_ventanas_total(i) = sum(idx);

    %% DF completa del paciente

    valores_comp = T.DF_completo_Hz(idx);
    valores_comp = valores_comp(isfinite(valores_comp));

    N_DF_completo(i) = numel(valores_comp);

    if ~isempty(valores_comp)

        DF_completo_media_Hz(i) = ...
            mean(valores_comp, 'omitnan');

        DF_completo_std_Hz(i) = ...
            std(valores_comp, 0, 'omitnan');
    end

    %% DF residual del paciente

    valores_res = T.DF_residual_Hz(idx);
    valores_res = valores_res(isfinite(valores_res));

    N_DF_residual(i) = numel(valores_res);

    if ~isempty(valores_res)

        DF_residual_media_Hz(i) = ...
            mean(valores_res, 'omitnan');

        DF_residual_std_Hz(i) = ...
            std(valores_res, 0, 'omitnan');
    end
end

T_pacientes = table( ...
    Grupo, ...
    Ritmo, ...
    Paciente, ...
    Base, ...
    Registro, ...
    N_ventanas_total, ...
    N_DF_completo, ...
    DF_completo_media_Hz, ...
    DF_completo_std_Hz, ...
    N_DF_residual, ...
    DF_residual_media_Hz, ...
    DF_residual_std_Hz);

writetable( ...
    T_pacientes, ...
    fullfile( ...
    carpeta_out, ...
    'resultados_DF_por_paciente.xlsx'));

%% ============================================================
% RESUMEN POR GRUPO A PARTIR DE LOS PACIENTES
%
% Cada paciente aporta una sola media.
%% ============================================================

resumen_grupos = table();

for g = 1:numel(grupos)

    grupo = string(grupos{g});

    idx_grupo = T_pacientes.Grupo == grupo;

    %% DF completa por paciente

    x_comp = ...
        T_pacientes.DF_completo_media_Hz(idx_grupo);

    x_comp = x_comp(isfinite(x_comp));

    %% DF residual por paciente

    x_res = ...
        T_pacientes.DF_residual_media_Hz(idx_grupo);

    x_res = x_res(isfinite(x_res));

    %% Crear fila resumen

    fila = table( ...
        grupo, ...
        sum(idx_grupo), ...
        numel(x_comp), ...
        string(sprintf( ...
        '%.3f +/- %.3f', ...
        mean(x_comp, 'omitnan'), ...
        std(x_comp, 0, 'omitnan'))), ...
        numel(x_res), ...
        string(sprintf( ...
        '%.3f +/- %.3f', ...
        mean(x_res, 'omitnan'), ...
        std(x_res, 0, 'omitnan'))), ...
        'VariableNames', { ...
        'Grupo', ...
        'N_pacientes_total', ...
        'N_pacientes_DF_completo', ...
        'DF_completo_Hz', ...
        'N_pacientes_DF_residual', ...
        'DF_residual_Hz'});

    resumen_grupos = ...
        [resumen_grupos; fila]; %#ok<AGROW>
end

writetable( ...
    resumen_grupos, ...
    fullfile( ...
    carpeta_out, ...
    'resumen_DF_por_grupo_pacientes.xlsx'));

disp(' ');
disp('RESUMEN POR GRUPO CALCULADO ENTRE PACIENTES');
disp(resumen_grupos);

%% ============================================================
% TABLAS FILTRADAS PARA GENERAR FIGURAS
%% ============================================================

T_comp_pacientes = T_pacientes( ...
    isfinite(T_pacientes.DF_completo_media_Hz), :);

T_res_pacientes = T_pacientes( ...
    isfinite(T_pacientes.DF_residual_media_Hz), :);

%% ============================================================
% BOXPLOTS POR PACIENTE
%% ============================================================

%% 1) DF completa: SANO frente a FA_PAROXISTICA_RS

idx_RS_comp = ...
    T_comp_pacientes.Grupo == "SANO" | ...
    T_comp_pacientes.Grupo == "FA_PAROXISTICA_RS";

if any(idx_RS_comp)

    f1 = figure( ...
        'Visible', 'off', ...
        'Color', 'w', ...
        'Position', [100 100 950 650]);

    boxplot( ...
        T_comp_pacientes.DF_completo_media_Hz(idx_RS_comp), ...
        T_comp_pacientes.Grupo(idx_RS_comp));

    title( ...
        'DF media del ECG completo por paciente en ritmo sinusal', ...
        'FontSize', 18, ...
        'FontWeight', 'bold');

    ylabel( ...
        'Frecuencia dominante media (Hz)', ...
        'FontSize', 20);

    xlabel('Grupo', ...
        'FontSize', 20);

    ylim([0.5 2]);

    set(gca, ...
        'FontSize', 18, ...
        'LineWidth', 1.2);

    grid on
    box on

    saveas( ...
        f1, ...
        fullfile( ...
        carpeta_out, ...
        'boxplot_DF_completo_por_paciente_RS.png'));

    close(f1);
end

%% 2) DF residual: SANO frente a FA_PAROXISTICA_RS

idx_RS_res = ...
    T_res_pacientes.Grupo == "SANO" | ...
    T_res_pacientes.Grupo == "FA_PAROXISTICA_RS";

if any(idx_RS_res)

    f2 = figure( ...
        'Visible', 'off', ...
        'Color', 'w', ...
        'Position', [100 100 950 650]);

    boxplot( ...
        T_res_pacientes.DF_residual_media_Hz(idx_RS_res), ...
        T_res_pacientes.Grupo(idx_RS_res));

    title( ...
        'DF media del ECG residual por paciente en ritmo sinusal', ...
        'FontSize', 18, ...
        'FontWeight', 'bold');

    ylabel( ...
        'Frecuencia dominante media (Hz)', ...
        'FontSize', 20);

    xlabel('Grupo', ...
        'FontSize', 20);

    ylim([0.5 2]);

    set(gca, ...
        'FontSize', 18, ...
        'LineWidth', 1.2);

    grid on
    box on

    saveas( ...
        f2, ...
        fullfile( ...
        carpeta_out, ...
        'boxplot_DF_residual_por_paciente_RS.png'));

    close(f2);
end

%% 3) DF completa: FA_PAROXISTICA_FA frente a FA_PERSISTENTE

idx_FA_comp = ...
    T_comp_pacientes.Grupo == "FA_PAROXISTICA_FA" | ...
    T_comp_pacientes.Grupo == "FA_PERSISTENTE";

if any(idx_FA_comp)

    f3 = figure( ...
        'Visible', 'off', ...
        'Color', 'w', ...
        'Position', [100 100 950 650]);

    boxplot( ...
        T_comp_pacientes.DF_completo_media_Hz(idx_FA_comp), ...
        T_comp_pacientes.Grupo(idx_FA_comp));

    title( ...
        'DF media del ECG completo por paciente en fibrilación auricular', ...
        'FontSize', 18, ...
        'FontWeight', 'bold');

    ylabel( ...
        'Frecuencia dominante media (Hz)', ...
        'FontSize', 20);

    xlabel('Grupo', ...
        'FontSize', 20);

    ylim([3 9]);

    set(gca, ...
        'FontSize', 18, ...
        'LineWidth', 1.2);

    grid on
    box on

    saveas( ...
        f3, ...
        fullfile( ...
        carpeta_out, ...
        'boxplot_DF_completo_por_paciente_FA.png'));

    close(f3);
end

%% 4) DF residual: FA_PAROXISTICA_FA frente a FA_PERSISTENTE

idx_FA_res = ...
    T_res_pacientes.Grupo == "FA_PAROXISTICA_FA" | ...
    T_res_pacientes.Grupo == "FA_PERSISTENTE";

if any(idx_FA_res)

    f4 = figure( ...
        'Visible', 'off', ...
        'Color', 'w', ...
        'Position', [100 100 950 650]);

    boxplot( ...
        T_res_pacientes.DF_residual_media_Hz(idx_FA_res), ...
        T_res_pacientes.Grupo(idx_FA_res));

    title( ...
        'DF media del ECG residual por paciente en fibrilación auricular', ...
        'FontSize', 18, ...
        'FontWeight', 'bold');

    ylabel( ...
        'Frecuencia dominante media (Hz)', ...
        'FontSize', 20);

    xlabel('Grupo', ...
        'FontSize', 20);

    ylim([3 9]);

    set(gca, ...
        'FontSize', 18, ...
        'LineWidth', 1.2);

    grid on
    box on

    saveas( ...
        f4, ...
        fullfile( ...
        carpeta_out, ...
        'boxplot_DF_residual_por_paciente_FA.png'));

    close(f4);
end

%% ============================================================
% HISTOGRAMAS POR GRUPO Y POR PACIENTE
%% ============================================================

bordes_RS = 0.5:0.05:2;
bordes_FA = 3:0.10:9;

for g = 1:numel(grupos)

    grupo = string(grupos{g});

    x_comp = ...
        T_pacientes.DF_completo_media_Hz( ...
        T_pacientes.Grupo == grupo);

    x_res = ...
        T_pacientes.DF_residual_media_Hz( ...
        T_pacientes.Grupo == grupo);

    x_comp = x_comp(isfinite(x_comp));
    x_res = x_res(isfinite(x_res));

    %% Elegir bordes según el ritmo

    if grupo == "SANO" || ...
            grupo == "FA_PAROXISTICA_RS"

        bordes = bordes_RS;
        limites_x = [0.5 2];

    elseif grupo == "FA_PAROXISTICA_FA" || ...
            grupo == "FA_PERSISTENTE"

        bordes = bordes_FA;
        limites_x = [3 9];

    else

        warning( ...
            'Grupo no reconocido para el histograma: %s', ...
            grupo);

        continue
    end

    %% Histograma DF completa por paciente

    if numel(x_comp) >= 2

        f = figure( ...
            'Visible', 'off', ...
            'Color', 'w', ...
            'Position', [100 100 950 650]);

        histogram(x_comp, bordes);

        title( ...
            ['DF media del ECG completo por paciente - ' ...
            char(grupo)], ...
            'Interpreter', 'none', ...
            'FontSize', 18, ...
            'FontWeight', 'bold');

        xlabel( ...
            'Frecuencia dominante media (Hz)', ...
            'FontSize', 20);

        ylabel( ...
            'Número de pacientes', ...
            'FontSize', 20);

        xlim(limites_x);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas( ...
            f, ...
            fullfile( ...
            carpeta_out, ...
            ['hist_DF_completo_por_paciente_' ...
            char(grupo) '.png']));

        close(f);
    end

    %% Histograma DF residual por paciente

    if numel(x_res) >= 2

        f = figure( ...
            'Visible', 'off', ...
            'Color', 'w', ...
            'Position', [100 100 950 650]);

        histogram(x_res, bordes);

        title( ...
            ['DF media del ECG residual por paciente - ' ...
            char(grupo)], ...
            'Interpreter', 'none', ...
            'FontSize', 18, ...
            'FontWeight', 'bold');

        xlabel( ...
            'Frecuencia dominante media (Hz)', ...
            'FontSize', 20);

        ylabel( ...
            'Número de pacientes', ...
            'FontSize', 20);

        xlim(limites_x);

        set(gca, ...
            'FontSize', 18, ...
            'LineWidth', 1.2);

        grid on
        box on

        saveas( ...
            f, ...
            fullfile( ...
            carpeta_out, ...
            ['hist_DF_residual_por_paciente_' ...
            char(grupo) '.png']));

        close(f);
    end
end

fprintf('\nFIN\n');
fprintf('Resultados guardados en:\n%s\n', carpeta_out);

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

[Pxx, f_axis] = pwelch( ...
    x, ...
    hamming(win), ...
    noverlap, ...
    Nfft, ...
    Fs);

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

%% Rechazar picos próximos a los bordes

if DF <= 3.1 || DF >= 8.9
    DF = NaN;
    peak_power = NaN;
    return
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

[Pxx, f_axis] = pwelch( ...
    x, ...
    hamming(win), ...
    noverlap, ...
    Nfft, ...
    Fs);

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

%% Rechazar picos próximos a los bordes

if DF <= 0.55 || DF >= 1.95
    DF = NaN;
    peak_power = NaN;
    return
end

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