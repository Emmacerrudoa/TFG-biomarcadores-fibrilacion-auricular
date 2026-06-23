clear
clc
close all

%% ============================================================
% WILCOXON PAREADO + CORRECCIÓN FDR
% ANÁLISIS TEMPORAL DE FA PAROXÍSTICA, VENTANAS DE 10 s
%
% Se procesan dos análisis:
%
%   1) Transiciones RS -> FA, ventanas de 10 s
%   2) Transiciones FA -> RS, ventanas de 10 s
%
% En cada análisis hay tres momentos:
%
%   MOMENTO 1: -180 a -170 s
%   MOMENTO 2:  -90 a  -80 s
%   MOMENTO 3:  -10 a    0 s
%
% Para cada análisis se realizan tres comparaciones pareadas:
%
%   Comparación A: Momento 2 - Momento 1
%   Comparación B: Momento 3 - Momento 1
%   Comparación C: Momento 3 - Momento 2
%
% En cada comparación se calcula:
%   - Media, DE, mediana e IQR en cada momento
%   - Diferencia: momento B - momento A
%   - Wilcoxon pareado mediante signrank
%   - Corrección FDR de Benjamini-Hochberg
%
% Cada ID_transicion debe aportar exactamente:
%   - Una fila con Orden_momento = 1
%   - Una fila con Orden_momento = 2
%   - Una fila con Orden_momento = 3
%
% La corrección FDR se aplica por separado en cada comparación.
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

analisis = { ...

    struct( ...
        'nombre', ...
        "TRANSICIONES RS -> FA, VENTANAS DE 10 s, 3 MOMENTOS", ...
        'archivo_entrada', ...
        'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_10s_RS_A_FA_3ventanas\biomarcadores_transiciones_RS_a_FA_10s_3ventanas.xlsx', ...
        'archivo_salida', ...
        'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_10s_RS_A_FA_3ventanas\resultados_Wilcoxon_FA_paroxistica_en_RS_10s_3ventanas_FDR.xlsx', ...
        'etiquetas_momentos', ...
        ["RS alejado (-180 a -170 s)", ...
         "RS intermedio (-90 a -80 s)", ...
         "RS previo (-10 a 0 s)"]), ...

    struct( ...
        'nombre', ...
        "TRANSICIONES FA -> RS, VENTANAS DE 10 s, 3 MOMENTOS", ...
        'archivo_entrada', ...
        'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_10s_FA_A_RS_3ventanas\biomarcadores_transiciones_FA_a_RS_10s_3ventanas.xlsx', ...
        'archivo_salida', ...
        'C:\Users\Emma\Documents\MATLAB\BIOMARCADORES_10s_FA_A_RS_3ventanas\resultados_Wilcoxon_FA_paroxistica_en_FA_10s_3ventanas_FDR.xlsx', ...
        'etiquetas_momentos', ...
        ["FA alejada (-180 a -170 s)", ...
         "FA intermedia (-90 a -80 s)", ...
         "FA previa (-10 a 0 s)"]) ...
};

%% ============================
% PROCESAR LOS DOS ANÁLISIS
%% ============================

for a = 1:numel(analisis)

    nombre = analisis{a}.nombre;
    archivo_entrada = analisis{a}.archivo_entrada;
    archivo_salida = analisis{a}.archivo_salida;
    etiquetas_momentos = analisis{a}.etiquetas_momentos;

    fprintf('\n====================================================\n');
    fprintf('%s\n', nombre);
    fprintf('====================================================\n');

    if ~isfile(archivo_entrada)

        warning('No se encuentra el archivo:\n%s', ...
            archivo_entrada);

        continue
    end

    analizar_transiciones_temporales_3ventanas( ...
        archivo_entrada, ...
        archivo_salida, ...
        nombre, ...
        etiquetas_momentos, ...
        alpha);
end

fprintf('\n====================================================\n');
fprintf('FIN DE TODOS LOS ANÁLISIS.\n');
fprintf('====================================================\n');

%% ============================================================
% FUNCIÓN PRINCIPAL
%% ============================================================

function analizar_transiciones_temporales_3ventanas( ...
    archivo_entrada, ...
    archivo_salida, ...
    nombre_analisis, ...
    etiquetas_momentos, ...
    alpha)

%% ============================
% LEER TABLA
%% ============================

T = readtable(archivo_entrada);

T.Properties.VariableNames = ...
    matlab.lang.makeValidName(T.Properties.VariableNames);

%% ============================
% COMPROBAR COLUMNAS OBLIGATORIAS
%% ============================

columnas_obligatorias = { ...
    'ID_transicion', ...
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

T.ID_transicion = ...
    strtrim(string(T.ID_transicion));

T.Momento = ...
    strtrim(string(T.Momento));

fprintf('\nArchivo cargado:\n%s\n', ...
    archivo_entrada);

fprintf('Filas totales: %d\n', ...
    height(T));

fprintf('Transiciones distintas: %d\n', ...
    numel(unique(T.ID_transicion)));

%% ============================
% COMPROBAR TRÍOS
%% ============================

ids = unique( ...
    T.ID_transicion, ...
    'stable');

ids_validos = strings(0,1);
ids_incompletos = strings(0,1);
ids_duplicados = strings(0,1);

for i = 1:numel(ids)

    id = ids(i);

    idx_id = ...
        T.ID_transicion == id;

    n_momento1 = sum( ...
        idx_id & ...
        T.Orden_momento == 1);

    n_momento2 = sum( ...
        idx_id & ...
        T.Orden_momento == 2);

    n_momento3 = sum( ...
        idx_id & ...
        T.Orden_momento == 3);

    if n_momento1 == 1 && ...
            n_momento2 == 1 && ...
            n_momento3 == 1

        ids_validos(end+1,1) = id; %#ok<SAGROW>

    elseif n_momento1 > 1 || ...
            n_momento2 > 1 || ...
            n_momento3 > 1

        ids_duplicados(end+1,1) = id; %#ok<SAGROW>

    else

        ids_incompletos(end+1,1) = id; %#ok<SAGROW>
    end
end

fprintf('\nControl de tríos:\n');

fprintf('  Transiciones completas con 3 momentos: %d\n', ...
    numel(ids_validos));

fprintf('  Transiciones incompletas: %d\n', ...
    numel(ids_incompletos));

fprintf('  Transiciones duplicadas: %d\n', ...
    numel(ids_duplicados));

if isempty(ids_validos)

    error(['No se encontraron transiciones con ' ...
        'los tres momentos completos en:\n%s'], ...
        archivo_entrada);
end

% Conservar únicamente transiciones completas y no duplicadas.
T = T( ...
    ismember(T.ID_transicion, ids_validos), ...
    :);

%% ============================
% IDENTIFICAR BIOMARCADORES
%% ============================

vars = T.Properties.VariableNames;

% Columnas de identificación, tiempo y control.
cols_excluir = { ...
    'Base', ...
    'Registro', ...
    'ID_transicion', ...
    'N_transicion_registro', ...
    'N_transicion_global', ...
    'Tiempo_inicio_RS_s', ...
    'Tiempo_inicio_FA_s', ...
    'Tiempo_fin_RS_s', ...
    'Tiempo_fin_FA_s', ...
    'Tiempo_transicion_RS_FA_s', ...
    'Tiempo_transicion_FA_RS_s', ...
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

    % Protección adicional para excluir variables de tiempo,
    % ventanas, archivos e identificadores.
    if startsWith(v, 'Tiempo_', 'IgnoreCase', true) || ...
            startsWith(v, 'Ventana_', 'IgnoreCase', true) || ...
            startsWith(v, 'Archivo_', 'IgnoreCase', true) || ...
            startsWith(v, 'N_transicion', 'IgnoreCase', true)

        continue
    end

    if isnumeric(T.(v))

        biomarcadores{end+1} = v; %#ok<SAGROW>
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
% COMPARACIONES ENTRE MOMENTOS
%% ============================

comparaciones = [ ...
    1 2; ...
    1 3; ...
    2 3];

nombres_comparaciones = [ ...
    "Momento 2 - Momento 1"; ...
    "Momento 3 - Momento 1"; ...
    "Momento 3 - Momento 2"];

Resultados_total = table;

for c = 1:size(comparaciones,1)

    momento_a = comparaciones(c,1);
    momento_b = comparaciones(c,2);

    etiqueta_a = etiquetas_momentos(momento_a);
    etiqueta_b = etiquetas_momentos(momento_b);

    nombre_comparacion = nombres_comparaciones(c);

    fprintf('\n----------------------------------------------------\n');
    fprintf('%s\n', nombre_comparacion);
    fprintf('%s vs %s\n', etiqueta_b, etiqueta_a);
    fprintf('----------------------------------------------------\n');

    Resultados = table;

    %% ============================
    % WILCOXON PAREADO
    %% ============================

    for b = 1:numel(biomarcadores)

        biom = biomarcadores{b};

        x_a = ...
            NaN(numel(ids_validos),1);

        x_b = ...
            NaN(numel(ids_validos),1);

        %% CONSTRUIR LAS PAREJAS

        for i = 1:numel(ids_validos)

            id = ids_validos(i);

            idx_a = ...
                T.ID_transicion == id & ...
                T.Orden_momento == momento_a;

            idx_b = ...
                T.ID_transicion == id & ...
                T.Orden_momento == momento_b;

            x_a(i) = ...
                T.(biom)(idx_a);

            x_b(i) = ...
                T.(biom)(idx_b);
        end

        %% ELIMINAR NaN O Inf DE FORMA PAREADA

        idx_finito = ...
            isfinite(x_a) & ...
            isfinite(x_b);

        x1 = ...
            x_a(idx_finito);

        x2 = ...
            x_b(idx_finito);

        % Diferencia: momento B - momento A.
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

            media1 = ...
                mean(x1, 'omitnan');

            de1 = ...
                std(x1, 0, 'omitnan');

            mediana1 = ...
                median(x1, 'omitnan');

            iqr1 = ...
                iqr(x1);

            media2 = ...
                mean(x2, 'omitnan');

            de2 = ...
                std(x2, 0, 'omitnan');

            mediana2 = ...
                median(x2, 'omitnan');

            iqr2 = ...
                iqr(x2);

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

                significativo_sin_corregir = "Sí";

            else

                significativo_sin_corregir = "No";
            end

            %% DIRECCIÓN DEL CAMBIO

            if mediana_delta > 0

                direccion = ...
                    etiqueta_b + " > " + etiqueta_a;

            elseif mediana_delta < 0

                direccion = ...
                    etiqueta_b + " < " + etiqueta_a;

            else

                direccion = ...
                    etiqueta_b + " = " + etiqueta_a;
            end
        end

        %% CREAR FILA

        nueva_fila = table;

        nueva_fila.Analisis = ...
            string(nombre_analisis);

        nueva_fila.Comparacion = ...
            string(nombre_comparacion);

        nueva_fila.Momento_A = ...
            momento_a;

        nueva_fila.Momento_B = ...
            momento_b;

        nueva_fila.Etiqueta_Momento_A = ...
            string(etiqueta_a);

        nueva_fila.Etiqueta_Momento_B = ...
            string(etiqueta_b);

        nueva_fila.Biomarcador = ...
            string(biom);

        nueva_fila.N_parejas = ...
            n_parejas;

        nueva_fila.Media_Momento_A = ...
            media1;

        nueva_fila.DE_Momento_A = ...
            de1;

        nueva_fila.Mediana_Momento_A = ...
            mediana1;

        nueva_fila.IQR_Momento_A = ...
            iqr1;

        nueva_fila.Media_Momento_B = ...
            media2;

        nueva_fila.DE_Momento_B = ...
            de2;

        nueva_fila.Mediana_Momento_B = ...
            mediana2;

        nueva_fila.IQR_Momento_B = ...
            iqr2;

        nueva_fila.Media_Delta_B_menos_A = ...
            media_delta;

        nueva_fila.DE_Delta_B_menos_A = ...
            de_delta;

        nueva_fila.Mediana_Delta_B_menos_A = ...
            mediana_delta;

        nueva_fila.IQR_Delta_B_menos_A = ...
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
    % ORDENAR RESULTADOS DE ESTA COMPARACIÓN
    %% ============================

    Resultados = sortrows( ...
        Resultados, ...
        {'p_ajustado_FDR','p_value'}, ...
        {'ascend','ascend'});

    %% ============================
    % RESUMEN EN CONSOLA
    %% ============================

    idx_nominal = ...
        Resultados.Significativo_sin_corregir == "Sí";

    idx_fdr = ...
        Resultados.Significativo_FDR == "Sí";

    fprintf('\nResumen de %s:\n', nombre_comparacion);

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
            'Mediana_Delta_B_menos_A', ...
            'p_value', ...
            'p_ajustado_FDR', ...
            'Direccion_Mediana', ...
            'Usar_UMAP_FDR'}));

    else

        fprintf('\nNo queda ningún biomarcador significativo tras FDR.\n');
    end

    Resultados_total = ...
        [Resultados_total; Resultados]; %#ok<AGROW>
end

%% ============================
% GUARDAR EXCEL
%% ============================

writetable( ...
    Resultados_total, ...
    archivo_salida);

fprintf('\nResultados guardados en:\n%s\n', ...
    archivo_salida);

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