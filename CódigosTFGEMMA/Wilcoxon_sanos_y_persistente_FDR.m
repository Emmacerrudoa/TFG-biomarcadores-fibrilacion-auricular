clear
clc
close all

%% ============================================================
% WILCOXON PAREADO + CORRECCIÓN FDR
% SUJETOS SANOS Y FA PERSISTENTE
%
% Comparación en ambos grupos:
%   Momento 1: minuto 1, de   0 a  60 s
%   Momento 2: minuto 3, de 120 a 180 s
%
% Cada paciente debe aportar exactamente:
%   - Una fila con Orden_momento = 1
%   - Una fila con Orden_momento = 2
%
% En cada grupo se calcula:
%   - Media, DE, mediana e IQR en cada momento
%   - Diferencia: minuto 3 - minuto 1
%   - Wilcoxon pareado mediante signrank
%   - Corrección FDR de Benjamini-Hochberg
%
% La corrección FDR se aplica por separado en cada grupo.
%
% Salidas:
%   resultados_Wilcoxon_sanos_60s_FDR.xlsx
%   resultados_Wilcoxon_FA_persistente_60s_FDR.xlsx
%
% IMPORTANTE:
%   - Las variables de identificación, tiempo y control no se analizan.
%   - La significación principal se interpreta mediante:
%
%         p_ajustado_FDR < 0.05
%
%   - Usar_UMAP_FDR marca los biomarcadores significativos tras FDR,
%     que después pueden utilizarse como punto de partida para UMAP.
%% ============================================================

%% ============================
% CONFIGURACIÓN
%% ============================

alpha = 0.05;

carpeta_sanos = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_60s_SANO';

archivo_sanos = fullfile( ...
    carpeta_sanos, ...
    'biomarcadores_sanos_60s.xlsx');

salida_sanos = fullfile( ...
    carpeta_sanos, ...
    'resultados_Wilcoxon_sanos_60s_FDR.xlsx');


carpeta_persistente = ...
    'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_60s_FA_PERSISTENTE';

archivo_persistente = fullfile( ...
    carpeta_persistente, ...
    'biomarcadores_FA_persistente_60s.xlsx');

salida_persistente = fullfile( ...
    carpeta_persistente, ...
    'resultados_Wilcoxon_FA_persistente_60s_FDR.xlsx');

%% ============================
% ANALIZAR SUJETOS SANOS
%% ============================

fprintf('\n====================================================\n');
fprintf('ANÁLISIS TEMPORAL EN SUJETOS SANOS\n');
fprintf('====================================================\n');

analizar_grupo_temporal( ...
    archivo_sanos, ...
    salida_sanos, ...
    "SANO", ...
    alpha);

%% ============================
% ANALIZAR FA PERSISTENTE
%% ============================

fprintf('\n====================================================\n');
fprintf('ANÁLISIS TEMPORAL EN FA PERSISTENTE\n');
fprintf('====================================================\n');

analizar_grupo_temporal( ...
    archivo_persistente, ...
    salida_persistente, ...
    "FA_PERSISTENTE", ...
    alpha);

fprintf('\n====================================================\n');
fprintf('FIN DE TODOS LOS ANÁLISIS.\n');
fprintf('====================================================\n');

%% ============================================================
% FUNCIÓN PRINCIPAL
%% ============================================================

function analizar_grupo_temporal( ...
    archivo_entrada, ...
    archivo_salida, ...
    nombre_grupo, ...
    alpha)

%% ============================
% COMPROBAR Y LEER TABLA
%% ============================

if ~isfile(archivo_entrada)

    error('No se encuentra el archivo:\n%s', ...
        archivo_entrada);
end

T = readtable(archivo_entrada);

T.Properties.VariableNames = ...
    matlab.lang.makeValidName(T.Properties.VariableNames);

%% ============================
% COMPROBAR COLUMNAS OBLIGATORIAS
%% ============================

columnas_obligatorias = { ...
    'Paciente', ...
    'Momento', ...
    'Orden_momento'};

for i = 1:numel(columnas_obligatorias)

    if ~ismember( ...
            columnas_obligatorias{i}, ...
            T.Properties.VariableNames)

        error('Falta la columna obligatoria %s en:\n%s', ...
            columnas_obligatorias{i}, ...
            archivo_entrada);
    end
end

%% ============================
% NORMALIZAR FORMATOS
%% ============================

T.Paciente = ...
    strtrim(string(T.Paciente));

T.Momento = ...
    strtrim(string(T.Momento));

fprintf('\nArchivo cargado:\n%s\n', ...
    archivo_entrada);

fprintf('Filas totales: %d\n', ...
    height(T));

fprintf('Pacientes distintos: %d\n', ...
    numel(unique(T.Paciente)));

%% ============================
% COMPROBAR PAREJAS
%% ============================

pacientes = unique( ...
    T.Paciente, ...
    'stable');

pacientes_validos = strings(0,1);
pacientes_incompletos = strings(0,1);
pacientes_duplicados = strings(0,1);

for i = 1:numel(pacientes)

    paciente = pacientes(i);

    idx_paciente = ...
        T.Paciente == paciente;

    n_momento1 = sum( ...
        idx_paciente & ...
        T.Orden_momento == 1);

    n_momento2 = sum( ...
        idx_paciente & ...
        T.Orden_momento == 2);

    if n_momento1 == 1 && ...
            n_momento2 == 1

        pacientes_validos(end+1,1) = ...
            paciente; %#ok<SAGROW>

    elseif n_momento1 > 1 || ...
            n_momento2 > 1

        pacientes_duplicados(end+1,1) = ...
            paciente; %#ok<SAGROW>

    else

        pacientes_incompletos(end+1,1) = ...
            paciente; %#ok<SAGROW>
    end
end

fprintf('\nControl de parejas:\n');

fprintf('  Pacientes completos: %d\n', ...
    numel(pacientes_validos));

fprintf('  Pacientes incompletos: %d\n', ...
    numel(pacientes_incompletos));

fprintf('  Pacientes duplicados: %d\n', ...
    numel(pacientes_duplicados));

if isempty(pacientes_validos)

    error(['No se encontraron pacientes con ' ...
        'los dos momentos completos en:\n%s'], ...
        archivo_entrada);
end

% Conservar únicamente pacientes completos y no duplicados.
T = T( ...
    ismember(T.Paciente, pacientes_validos), ...
    :);

%% ============================
% IDENTIFICAR BIOMARCADORES
%% ============================

vars = T.Properties.VariableNames;

% Columnas de identificación, tiempo y control.
cols_excluir = { ...
    'Base', ...
    'Registro', ...
    'Paciente', ...
    'Momento', ...
    'Orden_momento', ...
    'Ventana_ini_rel_s', ...
    'Ventana_fin_rel_s', ...
    'Tiempo_inicio_bloque_s', ...
    'Tiempo_fin_bloque_s', ...
    'Archivo_ventana_1', ...
    'Archivo_ventana_2', ...
    'N_R', ...
    'N_RR', ...
    'P_NumOndas', ...
    'T_NumOndas' ...
};

biomarcadores = {};

for i = 1:numel(vars)

    v = vars{i};

    if ismember(v, cols_excluir)
        continue
    end

    % Protección adicional para excluir automáticamente
    % variables de tiempo, ventanas, archivos e identificadores.
    if startsWith(v, 'Tiempo_', 'IgnoreCase', true) || ...
            startsWith(v, 'Ventana_', 'IgnoreCase', true) || ...
            startsWith(v, 'Archivo_', 'IgnoreCase', true) || ...
            startsWith(v, 'ID_', 'IgnoreCase', true)

        continue
    end

    if isnumeric(T.(v))

        biomarcadores{end+1} = ...
            v; %#ok<SAGROW>
    end
end

if isempty(biomarcadores)

    error('No se encontraron biomarcadores numéricos en:\n%s', ...
        archivo_entrada);
end

fprintf('\nBiomarcadores analizados: %d\n', ...
    numel(biomarcadores));

disp(string(biomarcadores(:)))

%% ============================
% WILCOXON PAREADO
%% ============================

Resultados = table;

for b = 1:numel(biomarcadores)

    biom = biomarcadores{b};

    x_minuto1 = ...
        NaN(numel(pacientes_validos),1);

    x_minuto3 = ...
        NaN(numel(pacientes_validos),1);

    %% CONSTRUIR LAS PAREJAS

    for i = 1:numel(pacientes_validos)

        paciente = pacientes_validos(i);

        idx_momento1 = ...
            T.Paciente == paciente & ...
            T.Orden_momento == 1;

        idx_momento2 = ...
            T.Paciente == paciente & ...
            T.Orden_momento == 2;

        x_minuto1(i) = ...
            T.(biom)(idx_momento1);

        x_minuto3(i) = ...
            T.(biom)(idx_momento2);
    end

    %% ELIMINAR NaN O Inf DE FORMA PAREADA

    idx_finito = ...
        isfinite(x_minuto1) & ...
        isfinite(x_minuto3);

    x1 = ...
        x_minuto1(idx_finito);

    x2 = ...
        x_minuto3(idx_finito);

    % Diferencia: minuto 3 - minuto 1.
    delta = ...
        x2 - x1;

    n_parejas = ...
        numel(delta);

    %% INICIALIZAR RESULTADOS

    media1 = NaN;
    de1 = NaN;
    mediana1 = NaN;
    iqr1 = NaN;

    media2 = NaN;
    de2 = NaN;
    mediana2 = NaN;
    iqr2 = NaN;

    media_delta = NaN;
    de_delta = NaN;
    mediana_delta = NaN;
    iqr_delta = NaN;

    p = NaN;

    significativo_sin_corregir = ...
        "No datos suficientes";

    direccion = "";

    %% ESTADÍSTICOS Y WILCOXON

    if n_parejas >= 3

        % Minuto 1.
        media1 = ...
            mean(x1, 'omitnan');

        de1 = ...
            std(x1, 0, 'omitnan');

        mediana1 = ...
            median(x1, 'omitnan');

        iqr1 = ...
            iqr(x1);

        % Minuto 3.
        media2 = ...
            mean(x2, 'omitnan');

        de2 = ...
            std(x2, 0, 'omitnan');

        mediana2 = ...
            median(x2, 'omitnan');

        iqr2 = ...
            iqr(x2);

        % Diferencia minuto 3 - minuto 1.
        media_delta = ...
            mean(delta, 'omitnan');

        de_delta = ...
            std(delta, 0, 'omitnan');

        mediana_delta = ...
            median(delta, 'omitnan');

        iqr_delta = ...
            iqr(delta);

        %% WILCOXON DE RANGOS CON SIGNO

        if any(delta ~= 0)

            p = ...
                signrank(x2, x1);

        else

            p = 1;
        end

        %% SIGNIFICACIÓN SIN CORREGIR

        if p < alpha

            significativo_sin_corregir = ...
                "Sí";

        else

            significativo_sin_corregir = ...
                "No";
        end

        %% DIRECCIÓN DEL CAMBIO

        if mediana_delta > 0

            direccion = ...
                "Minuto 3 > Minuto 1";

        elseif mediana_delta < 0

            direccion = ...
                "Minuto 3 < Minuto 1";

        else

            direccion = ...
                "Minuto 3 = Minuto 1";
        end
    end

    %% CREAR FILA

    nueva_fila = table;

    nueva_fila.Grupo = ...
        string(nombre_grupo);

    nueva_fila.Biomarcador = ...
        string(biom);

    nueva_fila.N_parejas = ...
        n_parejas;

    nueva_fila.Media_Minuto_1 = ...
        media1;

    nueva_fila.DE_Minuto_1 = ...
        de1;

    nueva_fila.Mediana_Minuto_1 = ...
        mediana1;

    nueva_fila.IQR_Minuto_1 = ...
        iqr1;

    nueva_fila.Media_Minuto_3 = ...
        media2;

    nueva_fila.DE_Minuto_3 = ...
        de2;

    nueva_fila.Mediana_Minuto_3 = ...
        mediana2;

    nueva_fila.IQR_Minuto_3 = ...
        iqr2;

    nueva_fila.Media_Delta = ...
        media_delta;

    nueva_fila.DE_Delta = ...
        de_delta;

    nueva_fila.Mediana_Delta = ...
        mediana_delta;

    nueva_fila.IQR_Delta = ...
        iqr_delta;

    nueva_fila.p_value = ...
        p;

    nueva_fila.Significativo_sin_corregir = ...
        string(significativo_sin_corregir);

    nueva_fila.Direccion_Mediana = ...
        string(direccion);

    Resultados = ...
        [Resultados; nueva_fila]; %#ok<AGROW>
end

%% ============================
% CORRECCIÓN FDR
% BENJAMINI-HOCHBERG
%% ============================

Resultados.p_ajustado_FDR = ...
    ajustar_p_benjamini_hochberg( ...
        Resultados.p_value);

Resultados.Significativo_FDR = ...
    strings(height(Resultados),1);

Resultados.Usar_UMAP_FDR = ...
    strings(height(Resultados),1);

for i = 1:height(Resultados)

    if isfinite( ...
            Resultados.p_ajustado_FDR(i))

        if Resultados.p_ajustado_FDR(i) < alpha

            Resultados.Significativo_FDR(i) = ...
                "Sí";

            Resultados.Usar_UMAP_FDR(i) = ...
                "Sí";

        else

            Resultados.Significativo_FDR(i) = ...
                "No";

            Resultados.Usar_UMAP_FDR(i) = ...
                "No";
        end

    else

        Resultados.Significativo_FDR(i) = ...
            "No datos suficientes";

        Resultados.Usar_UMAP_FDR(i) = ...
            "No";
    end
end

%% ============================
% ORDENAR
%% ============================

Resultados = sortrows( ...
    Resultados, ...
    {'p_ajustado_FDR','p_value'}, ...
    {'ascend','ascend'});

%% ============================
% GUARDAR EXCEL
%% ============================

writetable( ...
    Resultados, ...
    archivo_salida);

fprintf('\nResultados guardados en:\n%s\n', ...
    archivo_salida);

%% ============================
% RESUMEN EN CONSOLA
%% ============================

idx_nominal = ...
    Resultados.Significativo_sin_corregir == "Sí";

idx_fdr = ...
    Resultados.Significativo_FDR == "Sí";

fprintf('\nResumen de %s:\n', ...
    nombre_grupo);

fprintf('  Significativos sin corregir: %d de %d\n', ...
    sum(idx_nominal), ...
    height(Resultados));

fprintf('  Significativos tras FDR: %d de %d\n', ...
    sum(idx_fdr), ...
    height(Resultados));

if any(idx_fdr)

    fprintf('\nBiomarcadores significativos tras FDR:\n');

    disp(Resultados(idx_fdr, { ...
        'Biomarcador', ...
        'N_parejas', ...
        'Mediana_Delta', ...
        'p_value', ...
        'p_ajustado_FDR', ...
        'Direccion_Mediana', ...
        'Usar_UMAP_FDR'}));

else

    fprintf('\nNo queda ningún biomarcador significativo tras FDR.\n');
end

end

%% ============================================================
% FUNCIÓN BENJAMINI-HOCHBERG
%% ============================================================

function p_ajustado = ...
    ajustar_p_benjamini_hochberg(p)

p = ...
    double(p(:));

p_ajustado = ...
    NaN(size(p));

idx_validos = ...
    find(isfinite(p));

if isempty(idx_validos)
    return
end

p_validos = ...
    p(idx_validos);

[p_ordenados, orden] = ...
    sort(p_validos, 'ascend');

m = ...
    numel(p_ordenados);

q_ordenados = ...
    p_ordenados .* m ./ (1:m)';

% Garantizar monotonicidad desde el final.
for i = m-1:-1:1

    q_ordenados(i) = min( ...
        q_ordenados(i), ...
        q_ordenados(i+1));
end

% Los valores p ajustados no pueden superar 1.
q_ordenados( ...
    q_ordenados > 1) = 1;

% Recuperar el orden original.
q_validos = ...
    NaN(m,1);

q_validos(orden) = ...
    q_ordenados;

p_ajustado(idx_validos) = ...
    q_validos;

end