#**************************************************************************************#
#**************************************************************************************#
#
#                                 Escuela de Mortalidad                        
#                        Asociación Latinoamericana de Población
#                          Estimación de tablas de mortalidad
#
#         Fecha de creación:        26-11-2024
#         Actualizado por:          @HCastanheira, @APDataSc
#         Fecha de actualización:   08-12-2024
#         Institución:              Centro Latinoamericano y Caribeño de Demografía
#         Contacto:                 helena.cruz@cepal.org
#
#**************************************************************************************#
#**************************************************************************************#

# 0) Preámbulo ----

rm(list = ls())

library(data.table)
library(DemoToolsData)
library(dplyr)
library(u5mr)
library(DemoTools)
library(ggplot2)
library(MortalityLaws)
library(dplyr)

options(scipen = 999)
options(warn=-1)

source("scr/00_mort_functions.R")


"--------------------------------------------------------------------------------"

# A continuación se estiman las tablas de mortalidad para el 2022 a partir de 
# varias fuentes de información.  Una de las principales fuentes es el Censo 
# Ecuador 2022 que tiene como fecha de referencia el 30-11-2022.  
# Para más detalles: https://www.censoecuador.gob.ec/  


# 1) TM a partir de la mortalidad en la niñez de censos (U5MR) y tablas modelo ----

## 1.1 Carga de datos del Censo Ecuador - 2022 ----

(cmr_data <- data.table(
  agegrp = c(15, 20, 25, 30, 35, 40, 45), 
  women = c(725594, 728503, 671987, 623044, 591298, 558825, 476056),
  child_born = c(79100, 389089, 762811, 1069998, 1293875, 1393444, 1287071),
  child_dead = c(1757, 8282, 16627, 25003, 34509, 43514, 49648)
))


## 1.2 Estimación de q0_5 y q0_1 del censo 2022 - versión iussp ----
## https://demographicestimation.iussp.org/content/indirect-estimation-child-mortality

(ecu_nac22_qx_iussp <- u5mr_trussell_adj(
  cmr_data,
  women = "women",
  child_born = "child_born",
  child_dead = "child_dead",
  agegrp = "agegrp",
  model = "west",
  svy_year = 2022+10/12+30/365,
  sex = "both",
  variant = "iussp",
  e_0 = 70
))

write.csv(ecu_nac22_qx_iussp, file = "out/ecu_nac22_q0_5_iussp.csv") # exportando a la carpeta "out"


### U5MR - Gráfica
jpeg("out/UFMR_fig.jpg", width = 500, height = 350)

with(ecu_nac22_qx_iussp,
     plot(year, q5, type = "b", pch = 19,
          ylim = c(0.02, .034),
          col = "purple", xlab = "Fecha de referencia", ylab = "U5MR",
          main = paste0("Under-five mortality en Ecuador 2022 estimada usando\n",
                        "West CD y la versión de Trussell del método de Brass")))

with(ecu_nac22_qx_iussp, text(year, q5, agegrp, cex=0.65, pos=3,col="blue"))

legend("bottomleft", legend=c("U5MR Censo 2022"),
       col = c("purple"), lty = 1:1, cex=0.8)

dev.off()


## 1.3 Tabla de vida a partir de la mortalidad en la niñez ----
q05 <- ecu_nac22_qx_iussp$q5[2]*0.5 + ecu_nac22_qx_iussp$q5[3]*0.5 #promedio U5MR de 20-24 y 25-29  


### Estimación de tabla de vida de mujeres con un sólo parámetro U5MR y CD_West 
lt_Ec_oneP <- lt_model_cdun_match_single(type = "CD_West", 
                                 Sex = "f", 
                                 indicator = "5q0", 
                                 value = q05,
                                 OAnew = 100)
### Gráfica de q_x
(lt_OneP_fig <-  
lt_Ec_oneP %>% 
  ggplot() +
  geom_point(aes(Age, nqx)) +
  geom_line(aes(Age, nqx)) +
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,99)) +
  ggtitle("Probabilidad de muerte por edad simple a partir de la mortalidad en la \nniñez (5q0) del censo y tabla modelo CD familia West. Ecuador 2022") +
  theme_bw())

ggsave( plot = lt_OneP_fig,
        filename = 'out/lt_OneP_fig.jpg',
        height = 5,
        width  = 10 )

### Esperanza de vida al nacer para las mujeres
lt_Ec_oneP %>% setDT() %>% 
  .[Age==0, .(ex)] 

"---------------------------------FIN de sección---------------------------------"


# 2) TM a través de las defunciones del hogar de censos ----

## 2.1 Carga de los datos del módulo de mortalidad del Censo Ecuador 2022 ----
mort_hog <- fread("dat/CPV_Mortalidad_2022_Nacional.csv")
names(mort_hog) <- tolower(names(mort_hog)) # nombres de variables en minúscula


## 2.2 Preparación del numerador (defunciones) del módulo de mortalidad de hogares ----   

### Tabla de defunciones por edad y sexo
### Nota: se toman las defunciones ocurridas un año antes del levantamiento,
### es decir, desde noviembre de 2021 a ocubre de 2022.
(mort_hog <- mort_hog[  (m0202==2021 & m0201>=11) | (m0202==2022 & m0201<=10), 
                  .(deaths=.N), .(m0202, m0201, m04, m03)] %>%
                          .[ , 
                              `:=`( age = ifelse(m03>=100 & m03<999, 100, m03),
                                    m0201 = ifelse(m0201==99, 6, m0201) 
                                    ) ] %>% 
                          .[ , .(deaths = sum(deaths)), .(sex = m04, age)] %>% 
                              setorder(sex, age))

### Prorrateo de edades no declaradas (999)
(adj_mort_hog <- 
  mort_hog[ age!=999, 
           .( sex, age, deaths ) ] %>%
  .[ , p_deaths := deaths / sum( deaths ), .( sex ) ] %>%
  merge(
    mort_hog[ age==999,
             .( sex, na_deaths = deaths ) ],
    by = c( 'sex' )
  ) %>%
  .[ , deaths_adj := deaths + p_deaths * na_deaths ] %>%
  .[ , .( sex, age, deaths = deaths_adj ) ])

#### Comprobación
sum(mort_hog$deaths)
sum(adj_mort_hog$deaths)


### Gráfica de las defunciones de los hogares 2022
(deaths_hog_x1_fig <-
adj_mort_hog %>%
  ggplot( ) +
  geom_line( aes( x = age, y = deaths, color = factor( sex ), group = factor( sex ) ), 
             size = 1 ) +
  theme_bw())

ggsave( plot = deaths_hog_x1_fig,
        filename = 'out/deaths_hog_x1_fig.jpg',
        height = 5,
        width  = 10 )


### Dado el problema de mala declaración de edad de las defunciones se decide 
### trabajar por edad quinquenal
(deaths_hog <- adj_mort_hog[ , `:=`(age5 = ifelse(age==0, 0,
                                            ifelse(age>=1 & age<=4, 1, 
                                                   age - age %% 5)), 
                                    sex = ifelse(sex==1, "m", "f"))] %>% 
                .[ , .(deaths = sum(deaths)), .(sex, age5)])  

### Gráfica de defunciones por edad quinquenal
(deaths_hog_x5_fig <-  
  deaths_hog %>%
  ggplot( ) +
  geom_line( aes( x = age5, y = deaths, color = factor( sex ), group = factor( sex ) ), 
             size = 1 ) +
  theme_bw())

ggsave( plot = deaths_hog_x5_fig,
        filename = 'out/deaths_hog_x5_fig.jpg',
        height = 5,
        width  = 10 )


## 2.3 Preparación del denominador, población a mitad del período de las defunciones----   
## es decir, 30 de abril de 2022

### Carga de la población de los dos últimos censos 2010 y 2022
(pop_dt <- 
  fread( 'dat/ecu_pop_census1022.csv' ) %>%
  .[ , .(year, date_ref, sex, 
         age= ifelse(age==0, 0, ifelse(age %in% 1:4, 1, age - age %% 5)), 
         pop)] %>%
  setorder( date_ref, sex, age ))   


### Se acota la población a 100+ años y se agrega por edad quinquenal  
(pop_dt <- 
  pop_dt %>% copy %>%
  .[ , age := ifelse( age > 100, 100, age ) ] %>%
  .[ , 
     .( pop = sum( pop ) ),
     .( year, date_ref, sex, age ) ])

pop_dt[ , max(age),year] # 100+ - ok


### Se traslada la población al 30 de abril de 2022
pop_interp_1022 <- data.table()
for( s in c("m", "f") ){

  pop_interp_1022 <-
    rbind(
      pop_interp_1022,
      data.table(
            date_ref = 2022.5,
            sex = s,
            age5 = c(0, 1, seq(5, 100, 5)),
            pop_exp = exp_interp( t = 2022+3/12+30/365, 
                                 y1 = pop_dt[ sex == s & year == 2010 ]$pop, 
                                 y2 = pop_dt[ sex == s & year == 2022 ]$pop, 
                                 t2 = 2022+10/12+30/365, 
                                 t1 = 2010+10/12+28/365 )
      ))

}

pop_interp_1022[ , .(pop=sum(pop_exp)), .(sex)]


## 2.4 Preparacion input tabla de mortalidad: muertes hogares + IMR/U5MR + completitud + población----

### Completitud de los datos de defunciones - WPP 2022

corr_complt <- rbind(
  corr_complt <- 
    fread( 'dat/ECU_COMPLETITUD_DEFUNCIONES_WPP2022.csv' ),
  
  corr_complt[ date_ref==2021.5, .(date_ref = 2022.5, sex, vr_comp, vr_comp_new)]
)
corr_complt[ , vr_comp := vr_comp_new ]


### Consolidación input tabla mortalidad: muertes hogares + + población + IMR/U5MR + completitud  

q1_q5 <- fread("dat/ecu_q1_q5_1990_2050.csv")

(lt_input_2022h <- 
  merge(
    deaths_hog,
    pop_interp_1022,
    by = c( 'sex', 'age5' )) %>%
  merge(
    q1_q5[ , .( date_ref, sex, q1 = q0_1, q5 = q0_5 ) ],
    by = c( 'date_ref', 'sex' )
  ) %>%
  merge(
    corr_complt,
    by = c( 'date_ref', 'sex' )
  ) %>% 
  setorder( sex, age5 ))


## 2.5 Construción tabla de mortalidad x edad simple----
single_est <- TRUE # TRUE para edad simple o FALSE para grupos quinquenales

### Cálculo tasas de mortalidad corregidas por completitud
lt_input_2022h[ , mx := (deaths / pop_exp) / vr_comp_new]

### Tabla mortalidad nacional por edad simple para ambos sexos
lt_output_2022h <- data.table()
for( s in c( 'm', 'f' ) ){
  for( y in unique( lt_input_2022h$date_ref ) ){
    
    temp_dt <- lt_input_2022h[ sex == s & date_ref == y ]
    
    # utiliza datos de igme para ajustar las tasas de mortalidad < 1 y 1-4
    q0_1 <- unique( temp_dt$q1 ) / 1000
    q0_5 <- unique( temp_dt$q5 ) / 1000
    q1_4 <- 1 - ( 1 - q0_5 ) / ( 1 - q0_1 )
    
    new_a0_1 <- DemoTools::lt_rule_1a0_ak( q0 = q0_1, Sex = s )
    new_m0_1 <- q0_1 / ( 1 - ( 1 - new_a0_1 ) * q0_1 )
    new_a1_4 <- DemoTools::lt_rule_4a1_cd( M0 = new_m0_1, Sex = s, region = "w" )
    new_m1_4 <- q1_4 / ( 4 - ( 4 - new_a1_4 ) * q1_4 )
    temp_dt[ age5 == 0 ]$mx <- new_m0_1
    temp_dt[ age5 == 1 ]$mx <- new_m1_4
    
    temp_lt <-
      lt_ambiguous(
        nMx_or_nqx_or_lx = temp_dt$mx,
        type = "mx",
        Age = temp_dt$age5,
        Sex = s,
        a0rule = "cd", 
        region = "w",
        axmethod = "un",
        extrapLaw = 'Kannisto',
        extrapFrom = 75,
        extrapFit  = temp_dt[ age5 %in% 60:70 ]$age5,
        OAG = TRUE,
        OAnew = 100,
        radix = 1,
        Single = single_est ) %>%
      setDT %>%
      .[ , date_ref := y ] %>%
      .[ , sex := s ]
    
    lt_output_2022h <- 
      rbind(
        lt_output_2022h,
        temp_lt[ , .( lt_desc = 'TABMORT Mortalidad de Hogares del Censo 2022',
                      year = trunc( y ),
                      date_ref, 
                      sex,
                      age = Age, 
                      mx = round( nMx, 6 ), 
                      qx = round( nqx, 6 ),
                      ax = round( nAx, 6 ), 
                      lx = round( lx, 6 ), 
                      dx = round( ndx, 6 ), 
                      Lx = round( nLx, 6 ), 
                      Sx = round( Sx, 6 ),
                      Tx = round( Tx, 6 ), 
                      ex = round( ex, 6 ) ) ]
      )
    
  }
  
}

### Gráfica de qx
(qx_hog_fig <-
lt_output_2022h[ age < 100 ]  %>%
  ggplot( ) +
  geom_line( aes( x = age, y = qx, color = factor( sex ), group = factor(sex ) ), size = 1 ) +
  scale_y_log10()+
  labs(color='sex') +   
  theme_bw())

ggsave( plot = qx_hog_fig,
        filename = 'out/qx_hog_fig.jpg',
        height = 5,
        width  = 10 )

lt_output_2022h[ age == 0 ] %>% dcast( year ~ sex, value.var = 'ex' )

"---------------------------------FIN de sección---------------------------------"


# 3) TM a partir de los registros vitales y censos ----

### 3.1 Carga de las defunciones para el ano censal ----
(eevv_dt <- 
  fread( 'dat/eevv_dt.csv' ) %>%
  setorder( date_ref, sex, age ) %>%
  .[ year %in% c( 2022 ) ] %>%
  .[ , age := ifelse( age > 100, 100, age ) ] %>%
  .[ , 
     .( deaths = sum( deaths ) ),
     .( year, date_ref, sex, age ) ])

### 3.2 Preparación de los datos de defunciones para el ano censal ----

### Para cálculo del promedio para el año de referencia, si se usara más de un año
eevv_dt[ , year_new := ifelse( year %in% 2022, 2022, year ) ]

(mort_dt <- 
  eevv_dt[ , 
            .( deaths = mean( deaths ) ),
            .( year = year_new,
                    sex, age ) ] %>%
         .[ , .( year,
            date_ref = year + 0.5,
            sex, age, deaths ) ])

mort_dt[ , max( age ), year]
mort_dt[ ,.N, year ]


### 3.3 Poblacion a mitad de año del censo 2022 ----
pop_dt <- 
  fread( 'dat/ecu_pob_interp.csv' ) %>%
  .[ , .(date_ref = year + 0.5, sex, age, year, pop = pob)] %>% 
  setorder( date_ref, sex, age )  %>%
  .[ year %in% c( 2022 ) ] 

#### Edad maxima (grupo abierto) de los datos de poblacion - 100+ - ok
pop_dt[ , max( age, na.rm = TRUE ), year ]
pop_dt[ ,.N, year ]


### 3.4 Completitud de los datos de defunciones - wpp 2022 ----
corr_complt <- rbind(
  corr_complt <- 
    fread( 'dat/ECU_COMPLETITUD_DEFUNCIONES_WPP2022.csv' ),
  
  corr_complt[ date_ref==2021.5, .(date_ref = 2022.5, sex, vr_comp, vr_comp_new)]
)

corr_complt[ , vr_comp := vr_comp_new ]


### 3.5 Carga de la Mortalidad en la Ninez e Infantil de la Rev. 2024 del INEC ---- 
q1_q5 <- fread("dat/ecu_q1_q5_1990_2050.csv")


### 3.6 Preparacion input tabla mortalidad: muertes eevv + población + completitud + IMR/U5MR ----
lt_input <- 
  merge(
    mort_dt,
    pop_dt,
    by = c( 'year', 'date_ref', 'sex', 'age' )
  ) %>%
  merge(
    corr_complt,
    by = c( 'date_ref', 'sex' )
  ) %>%
  merge(
    q1_q5[ , .( date_ref, sex, q1 = q0_1, q5 = q0_5 ) ],
    by = c( 'date_ref', 'sex' )
  ) %>%
  setorder( date_ref, sex, age )


## 3.6 Cálculo y ajuste de las tasas de mortalidad----

### Cálculo de las tasas de mortalidad
lt_input[ , mx := deaths / pop  ]

### Suavización con Helligman & Pollard
lt_input[ ,
            mx_adj_hp := MortalityLaw(x = age,
                            mx = mx,
                            law = "HP", fit.this.x = 0:100,
                            opt.method = "LF2")$fitted.values,
  .( year, sex ) ]

### "Ajuste por completitud"
lt_input[ , mx_adj2 := mx_adj_hp / vr_comp]


### Gráfica de las m_x originales y ajustadas

(mx_eevv_censo_fig <-
lt_input %>%
  ggplot() + 
  geom_point(aes(x = age, y = mx, color = "mx"), size = 1) + 
  geom_line(aes(x = age, y = mx_adj_hp, color = "mx_adj_hp"), size = 1) + 
  geom_line(aes(x = age, y = mx_adj2, color = "mx_adj2"), size = 1, linetype = "dashed") + 
  scale_y_log10() + 
  scale_x_continuous(breaks = seq(0, 100, 10), labels = seq(0, 100, 10), limits = c(0, 100)) + 
  facet_grid(sex ~ year) + 
  theme_classic() +
  scale_color_manual(
    values = c("mx" = "orange", "mx_adj_hp" = "pink", "mx_adj2" = "blue"),
    name = "Tasas Ajustadas",
    breaks = c("mx", "mx_adj_hp", "mx_adj2"),
    labels = c("mx", "mx ajustada HP", "mx ajuste final")
  )) 


ggsave( plot = mx_eevv_censo_fig,
        filename = 'out/mx_eevv_censo_fig.jpg',
        height = 5,
        width  = 10 )


## 3.7 Construción tabla de mortalidad x edad simple----
single_est <- TRUE # TRUE para edad simples o FALSE para grupos quinquenales

lt_output <- data.table()
for( s in c( 'm', 'f' ) ){
  for( y in c( 2022.5 ) ){
    
    temp_dt <- lt_input[ sex == s & date_ref == y ]
    
    # utiliza datos de igme para ajustar las tasas de mortalidad < 1
    q0_1 <- unique( temp_dt$q1 ) / 1000
    new_a0_1 <- DemoTools::lt_rule_1a0_ak( q0 = q0_1, Sex = s )
    new_m0_1 <- q0_1 / ( 1 - ( 1 - new_a0_1 ) * q0_1 )
    temp_dt[ age == 0 ]$mx_adj2 <- new_m0_1
    
    temp_lt <- 
      lt_single_mx(
        nMx = temp_dt$mx_adj2,
        IMR = q0_1,
        Age = temp_dt$age,
        radix = 1,
        Sex = s,
        a0rule = 'ak',
        axmethod = "un",
        extrapLaw = 'Kannisto',
        extrapFrom = 75,
        extrapFit  = temp_dt[ age %in% 60:74 ]$age,
        OAG = TRUE,
        OAnew = 100
      ) %>%
      setDT %>%
      .[ , date_ref := y ] %>%
      .[ , sex := s ]
    
    lt_output <- 
      rbind(
        lt_output,
        temp_lt[ , .( lt_desc = 'TABMORT EEVV Censos (con ajuste)',
                      year = trunc( y ),
                      date_ref, 
                      sex,
                      age = Age, 
                      mx = round( nMx, 6 ), 
                      qx = round( nqx, 6 ),
                      ax = round( nAx, 6 ), 
                      lx = round( lx, 6 ), 
                      dx = round( ndx, 6 ), 
                      Lx = round( nLx, 6 ), 
                      Sx = round( Sx, 6 ),
                      Tx = round( Tx, 6 ), 
                      ex = round( ex, 6 ) ) ]
      )
    
  }
  
}

### Gráfica de qx
(qx_eevv_fig <-
    lt_output[ age < 100 ]  %>%
    ggplot( ) +
    geom_line( aes( x = age, y = qx, color = factor( sex ), group = factor(sex ) ), size = 1 ) +
    scale_y_log10()+
    labs(color='sex') +   
    theme_bw())

ggsave( plot = qx_eevv_fig,
        filename = 'out/qx_eevv_fig.jpg',
        height = 5,
        width  = 10 )


lt_output[ age == 0 ] %>% dcast( year ~ sex, value.var = 'ex' )



"---------------------------------FIN de sección---------------------------------"

# 4) TM a través de dos parámetros: mortalidad en la niñez (U5MR) y adulta (45q15) ----

## 4.1 Cálculo de q15_45 a partir de las EEVV y Censo (3) ----
(q15_45 <- 1 - lt_output[age==60 & sex=="f", .(lx)] %>% pull()/
              lt_output[age==15 & sex=="f", .(lx)] %>% pull())

## 4.2 Tabla de vida a partir de la mortalidad en la niñez y adulta----
lt_Ec_twoP <- lt_model_cdun_combin_single(type = "CD_West", 
                                                   Sex = "f", 
                                                   q1 = NA, 
                                                   q5 = q05, 
                                                   indicator_adult = "45q15", 
                                                   value_adult = q15_45, 
                                          OAnew = 100)

### Gráfica de q_x
(lt_TwoP_fig <-
lt_Ec_twoP %>% 
  ggplot() +
  geom_point(aes(Age, nqx)) +
  geom_line(aes(Age, nqx)) +
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,99)) +
  ggtitle("Probabilidad de muerte por edad simple a partir de la mortalidad en la niñez (5q0) de los censos \ny tablas modelo CD familia West. Ecuador 2022") +
  theme_bw())

ggsave( plot = lt_TwoP_fig,
        filename = 'out/lt_TwoP_fig.jpg',
        height = 5,
        width  = 10 )

### Esperanza de vida al nacer para las mujeres
lt_Ec_twoP %>% setDT() %>% 
  .[Age==0, .(ex)] 

"---------------------------------FIN de sección---------------------------------"


# 5) Análisis comparativo ----

lt_Ec_2022_WPP24 <- fread("dat/lt_ecu_2022_WPP24.csv")

## 5.1 Consolidación de estimaciones de distintas fuentes ----
lt <- rbind(
  lt_output_2022h[ sex=="f", .(source="Hogares", age, qx)],
  lt_output[ sex=="f", .(source="EEVV y Censo", age, qx)], 
  lt_Ec_oneP[ , .(source="U5MR Censo", age= Age, qx=nqx)],
  lt_Ec_twoP[ , .(source="U5MR Censo y 45q15", age= Age, qx=nqx)],
  lt_Ec_2022_WPP24[Sex=="Female", .(source="WPP 2024", age= AgeGrpStart, qx)]
  )

### Gráfica de q_x
(qx_Compara_fig <-
lt %>% 
  ggplot() +
  geom_line(data = . %>% filter(source == "U5MR Censo"), aes(age, qx, col = source)) +
  geom_point(data = . %>% filter(source == "EEVV y Censo"), aes(age, qx, col = source)) +
  geom_line(data = . %>% filter(source == "U5MR Censo y 45q15"), aes(age, qx, col = source)) +
  geom_line(data = . %>% filter(source == "Hogares"), aes(age, qx, col = source)) +
  geom_line(data = . %>% filter(source == "WPP 2024"), aes(age, qx, col = source)) +
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,99)) +
  ggtitle("Probabilidad de muerte de las mujeres por edad simple \na partir de las diversas fuentes. Ecuador 2022") +
  theme_bw())

ggsave( plot = qx_Compara_fig,
        filename = 'out/qx_Compara_fig.jpg',
        height = 5,
        width  = 10 )


## 5.2 Comparación de la e0 y U5MR de mujeres entre distintas fuentes ----

(resumen <- data.table(
Fuente=c("e0 con 1 parámetro (U5MR)", "e0 2 parámetros (U5MR y 45q15)", 
         "e0 hogares", "e0 censo y EEVV", "e0 WPP 2024"),
e0_2022 = round(c(
lt_Ec_oneP %>% setDT() %>% .[Age==0, .(ex)] %>% pull(),
lt_Ec_twoP %>% setDT() %>% .[Age==0, .(ex)] %>% pull(),
lt_output_2022h[ age == 0 & sex=="f" ] %>% dcast( year ~ sex, value.var = 'ex' ) %>% pull(),
lt_output[ age == 0 & sex=="f" ] %>% dcast( year ~ sex, value.var = 'ex' ) %>% pull(),
lt_Ec_2022_WPP24[ AgeGrpStart == 0 & Sex=="Female" ] %>% dcast( Time ~ Sex, value.var = 'ex' ) %>% pull()
), 2),
U5MRx1000 = round(c(1-lt_Ec_oneP[Age==5, .(lx)] %>% pull()/
                      lt_Ec_oneP[Age==0, .(lx)] %>% pull(),
         1-lt_Ec_twoP[Age==5, .(lx)] %>% pull()/
           lt_Ec_twoP[Age==0, .(lx)] %>% pull(),
         1-lt_output_2022h[age==5 & sex=="f", .(lx)] %>% pull()/
           lt_output_2022h[age==0 & sex=="f", .(lx)] %>% pull(),
         1-lt_output[age==5 & sex=="f", .(lx)] %>% pull()/
           lt_output[age==0 & sex=="f", .(lx)] %>% pull(),
         1-lt_Ec_2022_WPP24[AgeGrpStart==5 & Sex=="Female", .(lx)] %>% pull()/
           lt_Ec_2022_WPP24[AgeGrpStart==0 & Sex=="Female", .(lx)] %>% pull()
         )*1000, 2)
        )
       )

write.csv(resumen, file = "out/resumen.csv") # exportando a la carpeta "out"

"--------------------------------------------------------------------------------"
"--------------------------------------FIN---------------------------------------"
"--------------------------------------------------------------------------------"