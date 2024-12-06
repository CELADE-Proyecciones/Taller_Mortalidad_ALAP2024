#**************************************************************************************#
#**************************************************************************************#
#
#                                 Escuela de Mortalidad                        
#                        Asociación Latinoamericana de Población
#                          Estimación de tablas de mortalidad
#
#         Fecha de creación:        26-11-2024
#         Actualizado por:          @HCastanheira, @APDataSc
#         Fecha de actualización:   02-12-2024
#         Institución:              Centro Latinoamericano y Caribeño de Demografía
#         Contacto:                 helena.cruz@cepal.org
#
#**************************************************************************************#
#**************************************************************************************#

# 0) Preámbulo 

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

source("scr/00_mort_functions.R")


# A continuación se estiman las tablas de mortalidad para el 2022 a partir de varias fuentes 
# de información.  Una de las principales fuentes es el Censo Ecuador 2022 que tiene como 
# fecha de referencia el 30-11-2022.  Para más detalles: https://www.censoecuador.gob.ec/  


# 1) TM a partir de la mortalidad en la niñez de censos (U5MR) y tablas modelo ----

# Datos Censo Ecuador - 2022

cmr_data <- data.frame(
  agegrp = c(15, 20, 25, 30, 35, 40, 45), 
  women = c(725594, 728503, 671987, 623044, 591298, 558825, 476056),
  child_born = c(79100, 389089, 762811, 1069998, 1293875, 1393444, 1287071),
  child_dead = c(1757, 8282, 16627, 25003, 34509, 43514, 49648)
)


"q5 2022 - iussp"
ecu_nac22_qx_iussp <- u5mr_trussell_adj(
  cmr_data,
  women = "women",
  child_born = "child_born",
  child_dead = "child_dead",
  agegrp = "agegrp",
  model = "west",
  svy_year = 2022+10/12+30/365,
  sex = "both",
  variant = "iussp"
)

write.table(ecu_nac22_qx_iussp, "clipboard", row.names = F, sep = "\t")


## U5MR - Gráfica
with(ecu_nac22_qx_iussp,
     plot(year, q5, type = "b", pch = 19,
          ylim = c(0.02, .034),
          col = "purple", xlab = "Reference date", ylab = "u5MR",
          main = paste0("Under-five mortality, q(5) in Ecuador, estimated\n",
                        "using model West and the Trussell version of the Brass method")))

with(ecu_nac22_qx_iussp, text(year, q5, agegrp, cex=0.65, pos=3,col="blue"))


legend("bottomleft", legend=c("IUSSP"),
       col = c("purple"), lty = 1:1, cex=0.8)


# Tabla de vida a partir de la mortalidad en la niñez
q0_5 <- ecu_nac22_qx_iussp$q5[2]*0.5 + ecu_nac22_qx_iussp$q5[3]*0.5 
q01 <- ecu_nac22_qx_iussp$q1[2]*0.5 + ecu_nac22_qx_iussp$q1[3]*0.5 
year <- ecu_nac22_qx_iussp$year[2]*0.5 + ecu_nac22_qx_iussp$year[3]*0.5

lt_Ec_oneP <- lt_model_cdun_match_single(type = "CD_West", 
                                 Sex = "f", 
                                 indicator = "5q0", 
                                 value = q0_5,
                                 OAnew = 100)

lt_Ec_oneP %>% # Gráfica
  ggplot() +
  geom_point(aes(Age, nqx)) +
  geom_line(aes(Age, nqx)) +
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,100)) +
  ggtitle("Probabilidad de muerte por edad simple a partir de \nla mortalidad en la niñez (5q0) de los censos \ny tablas modelo CD familia West. Ecuador 2022") +
  theme_bw()


# Esperanza de vida al nacer para las mujeres
lt_Ec_oneP %>% setDT() %>% 
  .[Age==0, .(ex)] 

"---------------------------------FIN de sección---------------------------------"


# 2) TM a través de las defunciones del hogar de censos ----

mort_hog <- fread("dat/CPV_Mortalidad_2022_Nacional.csv")

names(mort_hog) <- tolower(names(mort_hog))

#Tabulado de edad y sexo
mort_hog <- mort_hog[  (m0202==2021 & m0201>=11) | (m0202==2022 & m0201<=10), 
                  .(deaths=.N), .(m0202, m0201, m04, m03)] %>%
                          .[ , 
                              `:=`( age = ifelse(m03>=100 & m03<999, 100, m03),
                                    m0201 = ifelse(m0201==99, 6, m0201) 
                                    ) ] %>% 
                          .[ , .(deaths = sum(deaths)), .(sex = m04, age)] %>% 
                              setorder(sex, age)

# Prorrateo de edades no declaradas
adj_mort_hog <- 
  mort_hog[ age!=999, 
           .( sex, age, deaths ) ] %>%
  .[ , p_deaths := deaths / sum( deaths ), .( sex ) ] %>%
  merge(
    mort_hog[ age==999,
             .( sex, na_deaths = deaths ) ],
    by = c( 'sex' )
  ) %>%
  .[ , deaths_adj := deaths + p_deaths * na_deaths ] %>%
  .[ , .( sex, age, deaths = deaths_adj ) ]

# Comprobación
sum(mort_hog$deaths)
sum(adj_mort_hog$deaths)

# Gráfica de las defunciones de los hogares 2022
adj_mort_hog %>%
  ggplot( ) +
  geom_line( aes( x = age, y = deaths, color = factor( sex ), group = factor( sex ) ), 
             size = 1 ) +
  theme_bw()

# Defunciones por edad quinquenal
deaths_hog <- adj_mort_hog[ , `:=`(age5 = ifelse(age==0, 0,
                                            ifelse(age>=1 & age<=4, 1, 
                                                   age - age %% 5)), 
                                    sex = ifelse(sex==1, "m", "f"))] %>% 
                .[ , .(deaths = sum(deaths)), .(sex, age5)]  

# Gráfica de defunciones por edad quinquenal
deaths_hog %>%
  ggplot( ) +
  geom_line( aes( x = age5, y = deaths, color = factor( sex ), group = factor( sex ) ), 
             size = 1 ) +
  theme_bw()



# Población de los dos últimos censos 2010 y 2022

pop_dt <- 
  fread( 'dat/ecu_pop_census1022.csv' ) %>%
  .[ , .(year, date_ref, sex, 
         age= ifelse(age==0, 0, ifelse(age %in% 1:4, 1, age - age %% 5)), 
         pop)] %>%
  setorder( date_ref, sex, age )   


# Ajusta poblacion > 100 
pop_dt <- 
  pop_dt %>% copy %>%
  .[ , age := ifelse( age > 100, 100, age ) ] %>%
  .[ , 
     .( pop = sum( pop ) ),
     .( year, date_ref, sex, age ) ]

pop_dt[ , max(age),year] # 100+ - ok



## Población al 30 de abril de 2022
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


# Tasas de mortalidad

# Preparacion input tabla mortalidad - muertes de hogares + población al 30 de abril 2022

q1_q5 <- fread("dat/ecu_q1_q5_1990_2050.csv")

lt_input_2022h <- 
  merge(
    deaths_hog,
    pop_interp_1022,
    by = c( 'sex', 'age5' )) %>%
  merge(
    q1_q5[ , .( date_ref, sex, q1 = q0_1, q5 = q0_5 ) ],
    by = c( 'date_ref', 'sex' )
  ) %>% 
  setorder( sex, age5 )


### Construcion tabla de mortalidad edad simple #---------------------------*

single_est <- TRUE # TRUE para edad simple o FALSE para grupos quinquenales


# Calcula tasas de mortalidad
lt_input_2022h[ , mx := deaths / pop_exp  ]


# Preparacion input tabla mortalidad nacional - eevv + censales 1974, 1982

lt_output_2022h <- data.table()
#tm_mx <- data.table()
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


lt_output_2022h[ age < 100 ]  %>%
  ggplot( ) +
  geom_line( aes( x = age, y = qx, color = factor( year ), group = factor(year ) ), size = 1 ) +
  scale_y_log10()+
  facet_wrap( ~ sex, ) +
  labs(color='año') +   
  # theme_bw()
  theme_classic()

lt_output_2022h[ age == 0 ] %>% dcast( year ~ sex, value.var = 'ex' )



"---------------------------------FIN de sección---------------------------------"


# 3) TM a partir de los registros vitales y censos ----

## 3.1 defunciones para el ano censal
eevv_dt <- 
  fread( 'dat/eevv_dt.csv' ) %>%
  setorder( date_ref, sex, age ) %>%
  .[ year %in% c( 2022 ) ] %>%
  .[ , age := ifelse( age > 100, 100, age ) ] %>%
  .[ , 
     .( deaths = sum( deaths ) ),
     .( year, date_ref, sex, age ) ]

# edad maxima (grupo abierto) de los datos de defunciones 95+ - ok
eevv_dt[ , max( age, na.rm = TRUE ), year ]

# distribucion de los missings - realizado por fuera de este script
adj_eevv_dt <- 
  eevv_dt

# calculo del promedio para el ano de referencia si se usara más de un año
adj_eevv_dt[ , year_new := ifelse( year %in% 2022, 2022, year ) ]

# datos preparados de defunciones
mort_dt <- 
  adj_eevv_dt[ , 
               .( deaths = mean( deaths ) ),
               .( year = year_new,
                  sex, age ) ] %>%
  .[ , .( year,
          date_ref = year + 0.5,
          sex, age, deaths ) ]

mort_dt[ , max( age ), year]
mort_dt[ ,.N, year ]


# 3.2 poblacion censada (en el medio del ano censal)
pop_dt <- 
  fread( 'dat/ecu_pob_interp.csv' ) %>%
  .[ , .(date_ref = year + 0.5, sex, age, year, pop = pob)] %>% 
  setorder( date_ref, sex, age )  %>%
  .[ year %in% c( 2022 ) ] 


pop_dt <- pop_dt[ , age := ifelse( age > 100, 100, age ) ] %>%
  .[ , 
     .( pop = sum( pop ) ),
     .( year, date_ref, sex, age ) ] 

# edad maxima (grupo abierto) de los datos de poblacion - 95+ - ok
pop_dt[ , max( age, na.rm = TRUE ), year ]
pop_dt[ ,.N, year ]

# 3.3 completitud de los datos de defunciones - wpp 2022

corr_complt <- rbind(
  corr_complt <- 
    fread( 'dat/ECU_COMPLETITUD_DEFUNCIONES_WPP2022.csv' ),
  
  corr_complt[ date_ref==2021.5, .(date_ref = 2022.5, sex, vr_comp, vr_comp_new)]
)

corr_complt[ , vr_comp := vr_comp_new ]


# Upload from file "ECU_Mort_Ninez_Infantil_Rev2024", sheet "Estimación U5MR e IMR" 
q1_q5 <- fread("dat/ecu_q1_q5_1990_2050.csv")


# 3.4 Preparacion input tabla mortalidad - eevv + censales 1990, 2001, 2010, 2019
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


# 3.5 Construción tabla de mortalidad edad simple #---------------------------*

single_est <- TRUE # TRUE para edad simples o FALSE para grupos quinquenales

# calcula tasas de mortalidad
lt_input[ , mx := deaths / pop  ]


# suavizacion con HP
lt_input[ ,
            mx_adj_hp := MortalityLaw(x = age,
                            mx = mx,
                            law = "HP", fit.this.x = 0:95,
                            opt.method = "LF2")$fitted.values,
  .( year, sex ) ]

# "Ajuste por completitud"
lt_input[ , mx_adj2 := mx_adj_hp / vr_comp]

"Gráfica"
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
  ) 


"Tablas de vida"
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


lt_output[ age == 0 ] %>% dcast( year ~ sex, value.var = 'ex' )




"---------------------------------FIN de sección---------------------------------"

# 4) TM a través de dos parámetros: mortalidad en la niñez y adulta (45q15) ----


q15_45 <- 1-lt_output[age==60 & sex=="f", .(lx)] %>% pull()/lt_output[age==15 & sex=="f", .(lx)] %>% pull()


lt_Ec_twoP <- lt_model_cdun_combin_single(type = "CD_West", 
                                                   Sex = "f", 
                                                   q1 = q01, 
                                                   q5 = NA, 
                                                   indicator_adult = "45q15", 
                                                   value_adult = q15_45, 
                                          OAnew = 100)

lt_Ec_twoP %>% # Gráfica
  ggplot() +
  geom_point(aes(Age, nqx)) +
  geom_line(aes(Age, nqx)) +
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,100)) +
  ggtitle("Probabilidad de muerte por edad simple a partir de \nla mortalidad en la niñez (5q0) de los censos \ny tablas modelo CD familia West. Ecuador 2022") +
  theme_bw()


lt_Ec_twoP %>% setDT() %>% 
  .[Age==0, .(ex)] 

"---------------------------------FIN de sección---------------------------------"


# 5) Análisis comparativos

lt <- rbind(
  lt_output_2022h[ sex=="f", .(source="Hogares", age, qx)],
  lt_output[ sex=="f", .(source="EEVV y Censo", age, qx)], 
  lt_Ec_oneP[ , .(source="U5MR Censo", age= Age, qx=nqx)],
  lt_Ec_twoP[ , .(source="U5MR Censo y 45q15", age= Age, qx=nqx)]
)


lt %>% 
  ggplot() +
  geom_point(data = . %>% filter(source == "U5MR Censo"), aes(age, qx, col = source)) +
  geom_point(data = . %>% filter(source == "EEVV y Censo"), aes(age, qx, col = source)) +
  geom_point(data = . %>% filter(source == "U5MR Censo y 45q15"), aes(age, qx, col = source)) +
  geom_point(data = . %>% filter(source == "Hogares"), aes(age, qx, col = source)) +
  scale_y_log10()+
  scale_x_continuous(breaks = seq(0,100,10), labels =  seq(0,100,10), 
                     limits = c(0,100)) +
  ggtitle("Probabilidad de muerte por edad simple a partir de \nla diversas fuentes. Ecuador 2022") +
  theme_bw()



# Comparación entre distintas fuentes

lt_Ec_oneP %>% setDT() %>% 
  .[Age==0, .(ex)] 

lt_output_2022h[ age == 0 & sex=="f" ] %>% dcast( year ~ sex, value.var = 'ex' )

lt_output[ age == 0 & sex=="f" ] %>% dcast( year ~ sex, value.var = 'ex' )

lt_Ec_twoP %>% setDT() %>% 
  .[Age==0, .(ex)] 
