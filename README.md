# Introducción a la evaluación, ajuste y análisis de datos de mortalidad en el contexto latinoamericano

<img src="img/logo_cepal.png" align="right" width="200" />

Repositorio del "Taller de estimación de la mortalidad en América Latina fuentes de datos y métodos" de la Escuela de Salud y Mortalidad de ALAP 2024 

### Centro Latinoamericano y Caribeño de Demografía (CELADE) - División de Población de la CEPAL
### Pre-evento del XI Congreso ALAP; Bogotá, Colombia
### 9 de diciembre de 2024 

### Facilitadora: Helena Cruz Castanheira
[helena.cruz@cepal.org](mailto:helena.cruz@cepal.org)

## Estructura de repositorio
* `dat` : datos de entrada
* `scr` : código para replicar el análisis y las visualizaciones
* `slides` : presentación del taller
* `out` : resultados


- [Primera parte: La estimación de la mortalidad en América Latina fuentes de datos y métodos](#primera-parte-la-estimación-de-la-mortalidad-en-américa-latina-fuentes-de-datos-y-métodos)
- [Segunda parte: Estimación de tablas de mortalidad con R a partir de varias fuentes](#Segunda-parte-Estimación-de-tablas-de-mortalidad-con-R-a-partir-de-varias-fuentes)
  - [TM a partir de la mortalidad en la niñez de censos (U5MR) y tablas modelo](#TM-a-partir-de-la-mortalidad-en-la-niñez-de-censos-U5MR-y-tablas-modelo)
  - [TM a través de las defunciones del hogar de censos](#TM-a-través-de-las-defunciones-del-hogar-de-censos)
  - [TM a partir de los registros vitales y censos](#TM-a-partir-de-los-registros-vitales-y-censos)
  - [TM a través de dos parámetros: mortalidad en la niñez y adulta (45q15)](#TM-a-través-de-dos-parámetros-mortalidad-en-la-niñez-y-adulta-45q15)
 
# Primera parte: La estimación de la mortalidad en América Latina fuentes de datos y métodos

Las diapositivas para la primera parte del taller están disponibles
[aquí](slides/ALAP_AnálisisMortalidad_v0.pdf).


# Segunda parte: Estimación de tablas de mortalidad con R a partir de varias fuentes

A continuación se estiman las tablas de mortalidad para el 2022 a partir de
varias fuentes de información.  Una de las principales fuentes es el Censo
Ecuador 2022 que tiene como fecha de referencia el 30-11-2022.
Para más detalles: https://www.censoecuador.gob.ec/

Antes de arrancar el taller práctico se recomienda tener instalado y cargado los siguientes paquetes:

``` r
library(data.table)
library(dplyr)
library(u5mr)
library(DemoTools)
library(DemoToolsData)
library(MortalityLaws)
library(ggplot2)
library(dplyr)
```

> Nota: si algún paquete no está instalado en su computadora, instálelo
> usando el comando `install.packages("dplyr")`.  En el caso del `DemoTools`:

``` r
# install.packages("remotes")

# requires the development version of rstan
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
remotes::install_github("timriffe/DemoTools")
```

## TM a partir de la mortalidad en la niñez de censos U5MR y tablas modelo

## TM a través de las defunciones del hogar de censos

## TM a partir de los registros vitales y censos

## TM a través de dos parámetros: mortalidad en la niñez y adulta 45q15
