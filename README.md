# TFG biomarcadores fibrilación auricular

Este repositorio contiene el código desarrollado para mi Trabajo de Fin de Grado, centrado en el análisis de señales ECG y la extracción de biomarcadores relacionados con la fibrilación auricular.

El objetivo principal del trabajo es estudiar diferencias entre sujetos sanos, pacientes con fibrilación auricular paroxística y pacientes con fibrilación auricular persistente, utilizando biomarcadores temporales, morfológicos y espectrales. Para ello, primero se calculan los biomarcadores, después se realiza el análisis estadístico y, finalmente, se emplea UMAP como herramienta de visualización para explorar los patrones observados.

## Contenido del repositorio

En este repositorio se incluyen scripts relacionados con:

* Lectura y preparación de registros ECG.
* Preprocesado y filtrado de la señal.
* Control de calidad de la señal y detección de artefactos.
* Detección y corrección de picos R.
* Segmentación de ventanas en ritmo sinusal y fibrilación auricular.
* Cancelación QRST para obtener el residual auricular.
* Cálculo de biomarcadores temporales, morfológicos y espectrales.
* Análisis estadístico de los biomarcadores obtenidos.
* Selección de biomarcadores no redundantes.
* Representación de los biomarcadores seleccionados mediante UMAP.
* Generación de tablas y figuras utilizadas en el trabajo.

## Bases de datos utilizadas

Para el desarrollo del trabajo se han utilizado registros ECG disponibles en PhysioNet:

* **MIT-BIH Normal Sinus Rhythm Database (NSRDB)** (`BASE_3_SANOS`).
* **MIT-BIH Atrial Fibrillation Database (AFDB)** (`BASE_2`).
* **Long-Term AF Database (LTAFDB)** (`BASE_1`).

Los datos originales no están incluidos en este repositorio. Para ejecutar el código es necesario descargar previamente las bases de datos desde PhysioNet y adaptar las rutas correspondientes en los scripts.

## Biomarcadores calculados

Los biomarcadores se agrupan en tres bloques principales.

### Biomarcadores temporales

A partir de la serie de intervalos RR se calculan diferentes parámetros de variabilidad del ritmo cardíaco, como:

* `RR_mean`
* `SDNN`
* `RMSSD`
* `SDSD`
* `pNN20`
* `pNN50`
* `CV_RR`
* `SD1`
* `SD2`
* `SD1_SD2_ratio`

### Biomarcadores morfológicos

En ritmo sinusal se analizan características de la onda P y de la onda T, como la amplitud, la dispersión y la correlación entre ondas de un mismo registro.

Algunos de los biomarcadores calculados son:

* `P_CorrIntraMedia`
* `P_CorrIntraStd`
* `P_AmpMedia`
* `P_AmpStd`
* `P_StdMedia`
* `T_CorrIntraMedia`
* `T_CorrIntraStd`
* `T_AmpMedia`
* `T_AmpStd`
* `T_StdMedia`

### Biomarcadores espectrales

También se calcula la frecuencia dominante de la señal ECG o del residual auricular, dependiendo del tipo de segmento analizado.

Algunos ejemplos son:

* `DF_completo_Hz`
* `DF_residual_Hz`

Para este análisis se utiliza el método de Welch, con ventanas de 30 s, solapamiento del 50 % y una NFFT definida como la siguiente potencia de dos superior a la longitud de la ventana.

## Análisis realizados

El código permite realizar diferentes análisis incluidos en el TFG.

### Comparación entre grupos

Se comparan biomarcadores entre sujetos sanos y pacientes con fibrilación auricular paroxística en ritmo sinusal, así como entre pacientes con fibrilación auricular paroxística durante FA y pacientes con FA persistente.

### Análisis temporal

Se estudian cambios entre ventanas de un mismo grupo o episodio, por ejemplo entre el minuto 1 y el minuto 3.

### Análisis de transiciones

También se analizan ventanas próximas y alejadas a los cambios de ritmo:

* De ritmo sinusal a fibrilación auricular: `RS → FA`.
* De fibrilación auricular a ritmo sinusal: `FA → RS`.

### Análisis estadístico

Una vez calculados los biomarcadores, se realizan las comparaciones estadísticas correspondientes. En función del análisis, se emplean pruebas pareadas o no pareadas, y se aplica corrección por comparaciones múltiples.

Los análisis estadísticos utilizados incluyen:

* Prueba de Wilcoxon para comparaciones pareadas.
* Prueba de Mann-Whitney para comparaciones entre grupos independientes.
* Corrección FDR de Benjamini-Hochberg para comparaciones múltiples.

### Selección de biomarcadores

Tras el análisis estadístico, se seleccionan biomarcadores representativos para evitar incluir variables muy redundantes entre sí. Esta selección se utiliza especialmente para las representaciones mediante UMAP.

### Visualización con UMAP

Finalmente, se utiliza UMAP para representar los biomarcadores seleccionados en un espacio de menor dimensión. Esta visualización permite explorar si existen patrones o agrupaciones entre los grupos estudiados o entre distintas ventanas temporales.

## Requisitos

El código se ha desarrollado principalmente en **MATLAB**.

Para ejecutar los scripts pueden ser necesarias las siguientes herramientas:

* MATLAB.
* Signal Processing Toolbox.
* Statistics and Machine Learning Toolbox.
* Registros ECG descargados desde PhysioNet.
* Código o funciones externas de UMAP añadidas al path de MATLAB, si se utilizan.

## Organización recomendada

Una posible organización del repositorio es:

## Organización recomendada

Una posible organización del repositorio, siguiendo el orden general del análisis realizado en el TFG, es:

```text
TFG-biomarcadores-fibrilación-auricular/
│
├── README.md
│
├── scripts/
│   ├── 01_preprocesado/
│   ├── 02_calidad_señal/
│   ├── 03_detección_R/
│   ├── 04_segmentación_ventanas/
│   ├── 05_cancelación_QRST/
│   ├── 06_cálculo_biomarcadores/
│   ├── 07_tablas_maestras/
│   ├── 08_análisis_estadístico/
│   ├── 09_selección_biomarcadores/
│   ├── 10_umap/
│   └── 11_visualización_figuras/
│
├── resultados/
│   ├── tablas/
│   ├── estadística/
│   ├── umap/
│   └── figuras/
│
└── documentación/
```

## Nota

Este repositorio se ha creado como apoyo al desarrollo del TFG y para organizar los scripts utilizados durante el análisis. Los registros ECG originales no se incluyen, por lo que deben descargarse directamente desde las bases de datos correspondientes.

## Autora

Emma Cerrudoa
Trabajo de Fin de Grado en Ingeniería Biomédica
