#**************************************************************************************#
#**************************************************************************************#
#
#                               Escuela de Mortalidad                        
#                        Asociación Latinoamericana de Población
#            Funciones auxiliares para la estimación de tablas de mortalidad
#
#         Fecha de creación:        26-11-2024
#         Actualizado por:          @HCastanheira, @APDataSc
#         Fecha de actualización:   02-12-2024
#         Institución:              Centro Latinoamericano y Caribeño de Demografía
#         Contacto:                 helena.cruz@cepal.org
#
#**************************************************************************************#
#**************************************************************************************#


# Model Life Table Implementations

# MATCH 
# --------------------------------------------------------------------------------------------------------------------------*
#      MATCH is derived from MORTPAK software package and customized for Abacus
#               UNITED NATION SOFTWARE PACKAGE FOR MORTALITY MEASUREMENT     
# Description of procedure
# This version is designed to use the same method as used in SLON and uses a look up table provided by the calling program
# Source of lookup table is MortCast MLTlookup
# Sara Hertog modified to maintain method, but simplify syntax for ccmppWPP 2021 revision
# and match to MortCast model life tables graduated to single year of age
# --------------------------------------------------------------------------------------------------------------------------*

#' Estimate UN or Coale Demeney family model life tables, matching on one input parameter
#' 
#' @details xxx
#' @param sex Choose the sex of the population. 
#' #' The following options are available: \itemize{
#'   \item{\code{"b"}} -- Both sex; 
#'   \item{\code{"f"}} -- Females;
#'   \item{\code{"m"}} -- Males.
#'   }
#' @param type Choose the family of model life tables
#'  #' The following options are available: \itemize{
#'   \item{\code{"CD_East"}} -- Coale-Demeny East; 
#'   \item{\code{"CD_North"}} -- Coale-Demeny North;
#'   \item{\code{"CD_South"}} -- Coale-Demeny South;
#'   \item{\code{"CD_West"}} -- Coale-Demeny West;
#'   \item{\code{"UN_Chilean"}} -- UN Chilean;
#'   \item{\code{"UN_Far_Eastern"}} -- UN Far Eastern;
#'   \item{\code{"UN_General"}} -- UN General;
#'   \item{\code{"UN_Latin_American"}} -- UN Latin American;
#'   \item{\code{"UN_South_Asian"}} -- UN South Asian.
#'   }
#' @param indicator character. Life table indicator to match on
#'  #' The following options are available: \itemize{
#'   \item{\code{"1q0"}} -- Probability of dying between birth and age 1; 
#'   \item{\code{"5q0"}} -- Probability of dying between birth and age 5; 
#'   \item{\code{"35q15"}} -- Probability of dying between age 15 and age 50;
#'   \item{\code{"45q15"}} -- Probability of dying between age 15 and age 60;
#'   \item{\code{"e0"}} -- Life expectancy at  birth;
#' @param values numeric. Values of the indicators to match on;
#' @inheritParams lt_abridged
#' @return data.frame. with two columns: age, giving abridged age groups, and mx
#' with model age specific mortality rates for abridged age groups
#' @importFrom stats uniroot MortCast DemoTools
#' @examples 
#' 
#' 
#'  lt <- lt_model_cdun_match_single(type = "CD_West", Sex = "f", indicator = "5q0", value = 0.150)
#'  lt <- lt_model_cdun_match_single(type = "CD_North", Sex = "m", indicator = "45q15", value = 0.450)


lt_model_cdun_match_single <- function(type,
                                       indicator, 
                                       value, 
                                       radix = 1e+05, 
                                       a0rule = "cd", 
                                       Sex = "m", 
                                       IMR = NA, 
                                       mod = TRUE, 
                                       SRB = 1.05, 
                                       OAnew = 130)   {
  
  # parse MLT lookup table from MortCast according to type and sex
  sexcode <- ifelse(Sex == "m", 1, ifelse(Sex == "f", 2, NA))
  MLTlookup <- MortCast::MLT1Ylookup
  MLTlookup <- as.data.frame(MLTlookup[MLTlookup$type == type & MLTlookup$sex == sexcode,])
  MLTlookup$index <- MLTlookup$e0
  
  region = "w"
  if (tolower(substr(type,1,2)) == tolower("CD")) {
    region <- tolower(substr(type,4,4))
  }
  
  if (indicator == "e0") {
    mlts       <- MLTlookup[, c("type","sex","index","age","mx")]
    # here we recompute e0 associated with model life tables to ensure that DemoTools functions return the input e0
    mlts$level <- NA
    for (level0 in unique(mlts$index)) {
      mlts$level[mlts$index == level0] <- DemoTools::lt_single_mx(nMx = mlts$mx[mlts$index == level0], 
                                                                  Age = mlts$age[mlts$index == level0], 
                                                                  Sex = Sex, a0rule = a0rule, region = region)$ex[1]
    }
  }
  
  if (indicator == "1q0") {
    # compute q1 levels for model life tables
    mlts       <- MLTlookup[MLTlookup$age == 1, c("type","sex","index","lx")]
    mlts$level <- 1-(mlts$lx / 100000)
    mlts       <- merge(MLTlookup, mlts[,c("type","sex","index","level")], by = c("type","sex","index"))
    
  }
  
  if (indicator == "5q0") {
    # compute q5 levels for model life tables
    mlts       <- MLTlookup[MLTlookup$age == 5, c("type","sex","index","lx")]
    mlts$level <- 1-(mlts$lx / 100000)
    mlts       <- merge(MLTlookup, mlts[,c("type","sex","index","level")], by = c("type","sex","index"))
    
  }
  
  if (indicator == "35q15") {
    # compute 35q15 levels for model life tables
    mlts       <- MLTlookup[MLTlookup$age %in% c(15,50), c("type","sex","index","age","lx")]
    mlts       <- reshape(mlts, direction = "wide", timevar = "age", 
                          idvar = c("type","sex","index"), sep = "_")
    mlts$level <- 1 - (mlts$lx_50 / mlts$lx_15)
    mlts       <- merge(MLTlookup, mlts[,c("type","sex","index","level")], by = c("type","sex","index"))
  }
  
  if (indicator == "45q15") {
    # compute 45q15 levels for model life tables
    mlts       <- MLTlookup[MLTlookup$age %in% c(15,60), c("type","sex","index","age","lx")]
    mlts       <- reshape(mlts, direction = "wide", timevar = "age", 
                          idvar = c("type","sex","index"), sep = "_")
    mlts$level <- 1 - (mlts$lx_60 / mlts$lx_15)
    mlts       <- merge(MLTlookup, mlts[,c("type","sex","index","level")], by = c("type","sex","index"))
  }
  
  # sort by level and age
  mlts         <- mlts[order(mlts$level, mlts$age),]
  
  # identify the model life tables with levels just below and above the value to match
  lvls   <- unique(mlts$level)
  iord   <- which.min(abs(value - lvls)) # identify closest level to value (could be higher or lower)
  lower  <- ifelse(lvls[iord] <= value, iord, iord-1)
  higher <- lower + 1
  ## PG: deal with out of bound levels
  if (lower==length(lvls)) {
    higher <- lower
    lower  <- lower-1
  }
  
  # parse the two matched mlts
  mlt_match <- mlts[mlts$level %in% c(lvls[c(lower,higher)]), c("level","age","mx")]
  mlt_match <- reshape(mlt_match, direction = "wide", timevar = "level", 
                       idvar = "age", sep = "_")
  
  # interpolate log(mx) between the two matched mlts according to the position of the input value
  # relative to the two matched levels (replicates Abacus approach)
  pct_val <- (value - lvls[lower]) / (lvls[higher]-lvls[lower])
  mx_hat  <- exp((1.0-pct_val)*log(mlt_match[,2])+pct_val*log(mlt_match[,3]))  # m(x,n)
  
  # compute the life table
  lt_out <- DemoTools::lt_single_mx(nMx = mx_hat, 
                                    Sex = Sex, 
                                    a0rule = a0rule,
                                    OAG = TRUE,
                                    OAnew = OAnew,
                                    radix = radix, 
                                    region = region, 
                                    IMR = IMR, 
                                    mod = mod, 
                                    SRB = SRB)
  
  
  return(lt_out)
  
}



"-------------------------------------------------------------------------------"
# COMBIN
# --------------------------------------------------------------------------------------------------------------------------*
#      COMBIN is derived from MORTPAK software package and customized for Abacus
#               UNITED NATION SOFTWARE PACKAGE FOR MORTALITY MEASUREMENT   
# Description of procedure
# Get MLT by matching on q(15,n) to get set of q's
# If Coale-Demeny models are used, use above q as a user defined model and use type = "user_defined" in BESTFT
# For first 2 age groups convert to l1 and l5.  
# If l1 is not available, use MLT and match on l5 to get l1.
# When both l1 and l5 are available, use BESTFT first with 1q0 and 5q20 then with 4q1 and 5q20. Average
# both results to get 5q5 and 5q10.  Note that 5q20 above is from MLT.
# If l5 is not available, use BESTFT with 1q0 and 5q20 to get 4q1, 5q5 and 5q10. 
# Note that early age groups come from second component fit.
# Sara Hertog modified to maintain method, but simplify syntax for ccmppWPP 2021 revision
# --------------------------------------------------------------------------------------------------------------------------*

#' Estimate UN or Coale Demeney family model life tables, matching multiple input parameters
#' 
#' @details xxx
#' @param Sex Choose the sex of the population. 
#' #' The following options are available: \itemize{
#'   \item{\code{"b"}} -- Both sex; 
#'   \item{\code{"f"}} -- Females;
#'   \item{\code{"m"}} -- Males.
#'   }
#' @param type Choose the family of model life tables
#'  #' The following options are available: \itemize{
#'   \item{\code{"CD_East"}} -- Coale-Demeny East; 
#'   \item{\code{"CD_North"}} -- Coale-Demeny North;
#'   \item{\code{"CD_South"}} -- Coale-Demeny South;
#'   \item{\code{"CD_West"}} -- Coale-Demeny West;
#'   \item{\code{"UN_Chilean"}} -- UN Chilean;
#'   \item{\code{"UN_Far_Eastern"}} -- UN Far Eastern;
#'   \item{\code{"UN_General"}} -- UN General;
#'   \item{\code{"UN_Latin_American"}} -- UN Latin American;
#'   \item{\code{"UN_South_Asian"}} -- UN South Asian.
#'   }
#' @param q1 numeric. Probability of dying between birth and age 1.
#' @param q5 numeric. Probability of dying between birth and age 5.
#' @param indicator_adult character. Adult mortality life table indicator to match on
#'  #' The following options are available: \itemize{
#'   \item{\code{"35q15"}} -- Probability of dying between age 15 and age 50;
#'   \item{\code{"45q15"}} -- Probability of dying between age 15 and age 60.
#' @param value_adult numeric. Value of the adult mortality indicator to match on;
#' @inheritParams lt_abridged
#' @return data.frame. with two columns: age, giving abridged age groups, and mx
#' with model age specific mortality rates for abridged age groups
#' @importFrom stats uniroot MortCast DemoTools
#' @examples 
#' 
#'
# lt <- lt_model_cdun_combin_single(type = "CD_West", Sex = "m",
#                                  q1 = NA, q5 = .05299,
#                                  indicator_adult = "45q15", value_adult = .210,
#                                  OAnew = 130)
#'



lt_model_cdun_combin_single <- function(type, 
                                        Sex, 
                                        q1 = NA, 
                                        q5 = NA, 
                                        indicator_adult, 
                                        value_adult,
                                        OAnew = 110,
                                        radix = 1e+05, 
                                        axmethod = "un", 
                                        a0rule = "cd",
                                        IMR = NA, 
                                        mod = TRUE, 
                                        SRB = 1.05, 
                                        extrapLaw = "kannisto", 
                                        extrapFrom = OAnew-5, 
                                        extrapFit = seq(OAnew-30, OAnew-5, 5))   {
  
  #############################*
  ##############################*
  ###############################*
  # New syntax for Abacus combin
  
  # indicator_adult can be "35q15", "45q15"
  
  # first match on adult mortality indicator via lt_model_cdun_match function;
  lt_temp_adult_single <- lt_model_cdun_match_single(type = type,
                                                     Sex  = Sex,
                                                     indicator = indicator_adult,
                                                     value = value_adult,
                                                     OAnew = 130)
  
  # compute abridged lt that corresponds to single
  lt_temp_adult_abr <- DemoTools::lt_single2abridged(lx  = lt_temp_adult_single$lx,
                                                     nLx = lt_temp_adult_single$nLx,
                                                     ex  = lt_temp_adult_single$ex)
  
  
  # extract model prob of dying bw ages 20 and 25
  q5_20 <- lt_temp_adult_abr$nqx[lt_temp_adult_abr$Age==20] 
  
  # if type is CD family, then use above output 5qx pattern as the user-defined model pattern
  if (tolower(type) %in% tolower(c("CD_West","CD_East","CD_North","CD_South"))) {
    users_model_pattern <- lt_temp_adult_abr$nqx[1:18]
    type_combin <- "user_defined"
  } else {
    users_model_pattern <- NA
    type_combin <- type
  }
  
  # If q5 is not available, use BESTFT with 1q0 and 5q20 to get 4q1, 5q5 and 5q10. 
  if (!is.na(q1) & is.na(q5)) {
    # bestft on q1 and 5q20 to get 4q1, 5q5, and 5q10
    lt_temp_child <- lt_model_un_bestft(type = type_combin,
                                        Sex  = Sex,
                                        age_start_abridged = c(0,20),
                                        qx_abridged = c(q1,q5_20),
                                        user_pattern = users_model_pattern,
                                        lt_compute = FALSE)
    lt_temp_child <- lt_temp_child[lt_temp_child$bestft_components == 2,]
    
    # splice input q1, best_ft q for ages 1-4, 5-9, 10-14, and adult match for older ages
    qx_out <- c(q1,lt_temp_child$nqx[2:4],lt_temp_adult_abr$nqx[5:nrow(lt_temp_adult_abr)])
    
  }
  
  # If q1 is not available, use MLT and match on q5 to get q1.
  if (is.na(q1) & !is.na(q5)) {
    # match on q5 to get q1
    lt_temp_child <- lt_model_cdun_match_single(type = type,
                                                Sex  = Sex,
                                                indicator = "5q0",
                                                value = q5)
    q1 <- lt_temp_child$nqx[1]
  }
  
  # When both q1 and q5 are available, use BESTFT first with 1q0 and 5q20 then with 4q1 and 5q20. Average
  # both results to get 5q5 and 5q10.  Note that 5q20 above is from MLT.
  if (!is.na(q1) & !is.na(q5)) {
    l1 <- 100000*(1-q1)
    l5 <- 100000*(1-q5)
    q1_4 <- 1-(l5/l1)
    # bestft first with 1q0 and 5q20
    lt_out1 <- lt_model_un_bestft(type = type_combin,
                                  Sex = Sex,
                                  age_start_abridged = c(0,20),
                                  qx_abridged = c(q1,q5_20),
                                  user_pattern = users_model_pattern,
                                  lt_compute = FALSE)
    qx_out1 <- lt_out1[lt_out1$bestft_components == 2, "nqx"]
    # bestft second with 4q1 and 5q20
    lt_out2 <- lt_model_un_bestft(type = type_combin,
                                  Sex = Sex,
                                  age_start_abridged = c(1,20),
                                  qx_abridged = c(q1_4,q5_20),
                                  user_pattern = users_model_pattern,
                                  lt_compute = FALSE)
    qx_out2 <- lt_out2[lt_out2$bestft_components == 2, "nqx"]
    
    qx_out <- c(q1, q1_4, c((qx_out1+qx_out2)/2)[3:4], lt_temp_adult_abr$nqx[5:nrow(lt_temp_adult_abr)])
  }
  
  region = "w"
  if (tolower(substr(type,1,2)) == tolower("CD")) {
    region <- tolower(substr(type,4,4))
  }
  
  lt_out_abr <- DemoTools::lt_abridged(nqx = qx_out,
                                       Age = lt_temp_adult_abr$Age, 
                                       Sex = Sex, 
                                       axmethod = axmethod, 
                                       a0rule = a0rule, 
                                       region = region,
                                       OAnew = 130,
                                       mod = mod,
                                       extrapFrom = max(lt_temp_adult_abr$Age) - 5)
  
  lt_out_sng <- DemoTools::lt_abridged2single(Age = lt_out_abr$Age,
                                              nMx = lt_out_abr$nMx,
                                              radix = radix,
                                              a0rule = a0rule, 
                                              Sex = Sex,
                                              region = region,
                                              IMR = q1, 
                                              mod = mod, 
                                              SRB = SRB, 
                                              OAG = TRUE,
                                              OAnew = OAnew,
                                              extrapLaw = extrapLaw,
                                              extrapFit = extrapFit,
                                              extrapFrom = extrapFrom)
  
  
  return(lt_out_sng)
  
}


"-------------------------------------------------------------------------------"
#BESTFT
# --------------------------------------------------------------------------------------------------------------------------*
#       BESTFT is derived from MORTPAK software package and customized for Abacus
#                  UNITED NATION SOFTWARE PACKAGE FOR MORTALITY MEASUREMENT    
# Sara Hertog modified to maintain method, but simplify syntax for ccmppWPP 2021 revision
# --------------------------------------------------------------------------------------------------------------------------*

#' Identify best fit life table from UN model life table families fitting up to three components
#' 
#' @details xxx
#' @param Sex Choose the sex of the population. 
#' #' The following options are available: \itemize{
#'   \item{\code{"f"}} -- Females;
#'   \item{\code{"m"}} -- Males.
#'   }
#' @param type Choose the family of model life tables
#'  #' The following options are available: \itemize{
#'   \item{\code{"UN_Chilean"}} -- UN Chilean;
#'   \item{\code{"UN_Far_Eastern"}} -- UN Far Eastern;
#'   \item{\code{"UN_General"}} -- UN General;
#'   \item{\code{"UN_Latin_American"}} -- UN Latin American;
#'   \item{\code{"UN_South_Asian"}} -- UN South Asian;
#'   \item{\code{"user_defined"}} -- user defined mortality pattern.
#'   }
#' @param age_start_abridged numeric. Vector of start ages of abridged age groups in input qx
#' @param qx_abridged numeric. Vector of probabilities of dying for abridged age groups to be fit
#' @param user_pattern numeric. Vector of qx to be used as model pattern when type=="user_defined"
#' @param lt_compute logical. Should output include all life table columns.  If FALSE, processing is faster and only nqx is output. 
#' @inheritParams lt_abridged
#' @return data.frame. Life table values for abridged age groups.
#' @importFrom stats uniroot MortCast DemoTools
#' @examples 


lt_model_un_bestft<-function(type, 
                             Sex, 
                             age_start_abridged, 
                             qx_abridged, 
                             user_pattern = NA,
                             lt_compute = TRUE,
                             OAnew = 110,
                             radix = 1e+05, 
                             axmethod = "un", 
                             a0rule = "ak", 
                             region = "w", 
                             IMR = NA, 
                             mod = TRUE, 
                             SRB = 1.05, 
                             extrapLaw = "kannisto", 
                             extrapFrom = OAnew-5, 
                             extrapFit = (OAnew-30):(OAnew-5)) {
  
  ## Do some checking/validating of input arguments
  
  Sex <- tolower(substr(Sex,1,1))
  sex_check <- Sex %in% c("m","f")
  
  if (!sex_check) 
    stop("The Sex argument must be either 'm' or 'f'.")
  
  type_check <- tolower(type) %in% tolower(c("user_defined", "UN_Latin_American", "UN_Far_Eastern",
                                             "UN_Chilean", "UN_South_Asian", "UN_General"))
  
  if (!type_check) 
    stop("The type argument must be one of 'user_defined', 'UN_Latin_American', 'UN_Far_Eastern',
                            'UN_Chilean', 'UN_South_Asian', 'UN_General'.")
  
  # abridged ages
  age <- c(0,1,seq(5,80,5))
  age_check <- all(age_start_abridged %in% age)
  
  if (!age_check) 
    stop("The age argument must be starting ages for abridged age groups 0,1,5,10,15....")
  
  qx_check <- all(!is.na(qx_abridged) & qx_abridged != 0 )
  
  if (!qx_check) 
    stop("The qx_abridged argument must be all non_zero values.")
  
  if (tolower(type) == tolower("user_defined")) {
    user_check <- length(user_pattern) == 18 # must be length 18 for this implementation
  }
  
  
  
  #  Empirical patterns for five UN life table families
  #  males at places 1:18; females at places 19:36
  
  if (tolower(type) == tolower("UN_Latin_American")) {
    # UN-Latin American
    EMP  <- c(-1.12977,-1.49127,-2.13005,-2.40748,-2.21892,-2.01157, 
              -1.93591,-1.86961,-1.76133,-1.64220,-1.49651,-1.34160,-1.15720, 
              -.96945,-.74708,-.52259,-.29449,-.04031,-1.22452,-1.45667, 
              -2.13881,-2.46676,-2.31810,-2.14505,-2.03883,-1.93924,-1.83147, 
              -1.74288,-1.62385,-1.47924,-1.28721,-1.07443,-.83152,-.59239, 
              -.35970,-.08623)
  }
  if (tolower(type) == tolower("UN_Chilean")) {
    # UN-Chilean
    EMP  <- c(-1.04722,-1.81992,-2.42430,-2.52487,-2.24491,-2.02821, 
              -1.90923,-1.78646,-1.66679,-1.52497,-1.37807,-1.21929,-1.03819, 
              -.84156,-.63201,-.42070,-.21110,+.01163,-1.12557,-1.82378, 
              -2.52319,-2.63933,-2.38847,-2.20417,-2.09701,-1.99128,-1.87930, 
              -1.75744,-1.61558,-1.45886,-1.26115,-1.05224,-.80346,-.58202, 
              -.35093,-.10587)
  }
  if (tolower(type) == tolower("UN_South_Asian")) {
    # UN-South Asian
    EMP  <- c(-.97864,-1.24228,-2.01695,-2.44280,-2.35424,-2.27012, 
              -2.16833,-2.05942,-1.90053,-1.71213,-1.51120,-1.28493,-1.08192, 
              -.84671,-.62964,-.40229,-.19622,-.00129,-0.97055,-1.15424, 
              -1.93962,-2.36857,-2.19082,-2.09358,-2.04788,-1.95922,-1.87311, 
              -1.76095,-1.61425,-1.39012,-1.15515,-0.90816,-.68011,-.43231, 
              -.17489,0.05948)
  }
  if (tolower(type) == tolower("UN_Far_Eastern")) {
    # UN-Far Eastern
    EMP  <- c(-1.53473,-2.15035,-2.61442,-2.66392,-2.42326,-2.23095, 
              -2.15279,-2.05765,-1.89129,-1.68244,-1.47626,-1.23020,-1.02801, 
              -.77148,-.54696,-.32996,-.11911,0.10572,-1.42596,-1.95200, 
              -2.55653,-2.68018,-2.33095,-2.15952,-2.03377,-1.94554,-1.82299, 
              -1.69084,-1.52189,-1.33505,-1.13791,-0.93765,-.72718,-.50916, 
              -.28389,-.01285)
  }
  if (tolower(type) == tolower("UN_General")) {
    # UN-General
    EMP  <- c(-1.27638,-1.78957,-2.35607,-2.55527,-2.34263,-2.16193, 
              -2.09109,-2.00215,-1.86781,-1.70806,-1.52834,-1.33100,-1.12934, 
              -.91064,-.68454,-.45685,-.23002,0.00844,-1.35963,-1.77385, 
              -2.39574,-2.64549,-2.44766,-2.28991,-2.18850,-2.08535,-1.97231, 
              -1.84731,-1.69291,-1.50842,-1.30344,-1.08323,-.84402,-.59485, 
              -.34158,-.06493)
  }
  
  if (tolower(type) == tolower("user_defined")) {
    # bring in user defined pattern
    EMP <- 0.50*log(user_pattern/(1.0-user_pattern))
    EMP <- rep(EMP,2)
  }
  
  if (Sex == "m") {
    
    EMP <- EMP[1:18]
    
    VEC <- cbind(c(.23686,.36077,.33445,.30540,.28931,.28678,.27950,.28023,.26073,
                   .23626,.20794,.17804,.15136,.13217,.12243,.11457,.10445,.08878),
                 c(-.46007,-.68813,.06414,.12479,.24384,.10713,.06507,.03339,.02833,
                   .06473,.08705,.10620,.11305,.09467,.10809,.14738,.21037,.30918),
                 c(.09331,-.29269,-.47139,-.17403,.10715,.28842,.33620,.33692,.21354,
                   .15269,.06569,.00045,-.03731,-.10636,-.11214,-.22258,-.19631,-.38123))
  }
  if (Sex == "f") {
    
    EMP <- EMP[19:36]
    
    VEC <- cbind(c(.18289,.31406,.31716,.30941,.32317,.32626,.30801,.29047,.25933,
                   .22187,.19241,.17244,.15729,.14282,.12711,.11815,.11591,.09772),
                 c(-.51009,-.52241,.08947,.03525,.03132,.07843,.06762,.00482,-.01409,
                   -.02178,.01870,.04427,.08201,.08061,.15756,.24236,.30138,.50530),
                 c(.23944,-.11117,.07566,.06268,-.26708,-.39053,-.28237,-.14277,-.05923,
                   .18909,.24773,.33679,.34121,.38290,.26731,.14442,.09697,-.13377))
  }
  
  
  
  # some parameters
  GAL1 <- VEC[,1]*VEC[,1]
  GAL2 <- VEC[,2]*VEC[,2]
  GAL3 <- VEC[,3]*VEC[,3]
  BAL1 <- VEC[,1]*VEC[,2]
  BAL2 <- VEC[,1]*VEC[,3]
  BAL3 <- VEC[,2]*VEC[,3]
  
  # initialize coefficients    
  GAMMA1 <- 0.0
  GAMMA2 <- 0.0
  GAMMA3 <- 0.0
  ALPHA1 <- 0.0
  ALPHA2 <- 0.0
  ALPHA3 <- 0.0
  BETA1 <- 0.0
  BETA2 <- 0.0
  BETA3 <- 0.0
  
  NCTR <- 0 # counter
  for (i in 1:length(age_start_abridged)){
    
    NCTR <- NCTR+1 # counter
    R      <- 0.5*log(qx_abridged[i]/(1.0-qx_abridged[i]))
    S      <- R-EMP[age == age_start_abridged[i]]
    SVEC1  <- S*VEC[age == age_start_abridged[i],1]
    SVEC2  <- S*VEC[age == age_start_abridged[i],2]
    SVEC3  <- S*VEC[age == age_start_abridged[i],3]
    GAMMA1 <- GAMMA1+GAL1[age == age_start_abridged[i]]
    GAMMA2 <- GAMMA2+GAL2[age == age_start_abridged[i]]
    GAMMA3 <- GAMMA3+GAL3[age == age_start_abridged[i]]
    ALPHA1 <- ALPHA1+SVEC1
    ALPHA2 <- ALPHA2+SVEC2
    ALPHA3 <- ALPHA3+SVEC3
    BETA1  <- BETA1+BAL1[age == age_start_abridged[i]]
    BETA2  <- BETA2+BAL2[age == age_start_abridged[i]]
    BETA3  <- BETA3+BAL3[age == age_start_abridged[i]]
    
  }
  
  lt_out <- list() # initialize component fits
  
  if(NCTR >= 3) NCTR <- 3
  NTEMP <- 4-NCTR
  for (i in NTEMP:3) {
    IC <- 4-i
    if(IC == 1 | IC == 2) GAMMA3 <- 1.0
    if(IC == 1 | IC == 2) BETA2 <- 0.0
    if(IC == 1 | IC == 2) BETA3 <- 0.0
    if(IC == 1) BETA1 <- 0.0
    if(IC == 1) GAMMA2 <- 1.0
    D  <- GAMMA1*GAMMA2*GAMMA3-GAMMA3*BETA1*BETA1-GAMMA2*BETA2*BETA2-GAMMA1*BETA3*BETA3+2.0*BETA1*BETA2*BETA3
    A1 <- ALPHA1*(GAMMA2*GAMMA3-BETA3*BETA3)+ALPHA2*(BETA2*BETA3-BETA1*GAMMA3)+ALPHA3*(BETA1*BETA3-BETA2*GAMMA2)
    A1 <- A1/D
    A2 <- ALPHA1*(BETA2*BETA3-BETA1*GAMMA3)+ALPHA2*(GAMMA1*GAMMA3-BETA2*BETA2)+ALPHA3*(BETA1*BETA2-BETA3*GAMMA1)
    A2 <- A2/D
    A3 <- ALPHA1*(BETA1*BETA3-BETA2*GAMMA2)+ALPHA2*(BETA1*BETA2-BETA3*GAMMA1)+ALPHA3*(GAMMA1*GAMMA2-BETA1*BETA1)
    A3 <- A3/D
    if (IC == 1) {
      A2 <- 0.0
      A3 <- 0.0
    }
    if (IC == 2) {
      A3 <- 0.0
    }
    CF <- EMP + A1*VEC[,1]+A2*VEC[,2]+A3*VEC[,3]
    CF <- exp(2.0*CF)/(1.0+exp(2.0*CF)) # nqx
    if (lt_compute == TRUE) {
      lt <- DemoTools::lt_abridged(Age = age, 
                                   nqx = CF, 
                                   Sex = Sex, 
                                   OAnew = OAnew,
                                   radix = radix, 
                                   axmethod = axmethod, 
                                   a0rule = a0rule, 
                                   region = region, 
                                   IMR = IMR, 
                                   mod = mod, 
                                   SRB = SRB, 
                                   extrapLaw = extrapLaw, 
                                   extrapFrom = extrapFrom, 
                                   extrapFit = extrapFit)
    } else {
      lt <- data.frame(Age = age, nqx = CF)
    }
    lt$bestft_components <- IC
    lt_out[[i]] <- lt
  }
  
  out.data <- do.call(rbind,lt_out)
  
  return(out.data) 
  
}  


# U5MR IUSSP & UN version function
u5mr_trussell_adj <- function (data, women = "women", child_born = "child_born", 
                               child_dead = "child_dead", agegrp = "agegrp", 
                               model = "west", svy_year = 1976.5, sex, 
                               variant = "iussp", e_0=60
) 
{
  agegrp <- data[[agegrp]]
  women <- data[[women]]
  child_born <- data[[child_born]]
  child_dead <- data[[child_dead]]
  pi <- child_born/women
  di <- child_dead/child_born
  p1 <- pi[1]
  p2 <- pi[2]
  p3 <- pi[3]
  coeff_trussell_ki <- coeff_trussell_ki
  t47 <- coeff_trussell_ki[coeff_trussell_ki$model == model, 
                           -1]
  ki <- t47$ai + (t47$bi * p1/p2) + (t47$ci * p2/p3)
  qx <- ki * di
  coeff_trussell_ti <- coeff_trussell_ti
  t48 <- coeff_trussell_ti[coeff_trussell_ti$model == model, 
                           -1]
  ti <- t48$ai + (t48$bi * p1/p2) + (t48$ci * p2/p3)
  year <- round(svy_year - ti, 1)
  
  if(variant=="iussp"){
    
    Ys <- DemoToolsData::modelLTx1 %>%
      setDT() %>%  
      .[type_mlt == "CD East" & e0 == e_0] %>% 
      dcast(age ~ sex, value.var = "lx1", fun.aggregate = sum) %>% 
      .[ , lx := (1.05*male+female)/2.05] %>% 
      .[ , Ys := 0.5*log((100000-lx)/lx)] %>% 
      .[ age %in% c(1:3, 5, 10, 15, 20) ] %>% 
      pull(Ys)
    
    alpha <- 0.5*(log(qx/(1-qx)))-Ys  
    
    q1 <- exp(2*(alpha+Ys[1]))/(1+(exp(2*(alpha+Ys[1]))))    
    q5 <- exp(2*(alpha+Ys[4]))/(1+(exp(2*(alpha+Ys[4]))))
    
  } else if (variant=="un") {
    
    coale_demeny_ltm <- coale_demeny_ltm
    inter <- coale_demeny_ltm[[paste0(model, "_", sex)]]
    age_index <- paste0("q", c(1, 2, 3, 5, 10, 15, 20))
    h <- vector("numeric", length = 7L)
    
    q5 <- sapply(1:7, function(x) {
      q <- qx[x]
      c <- inter[[age_index[x]]]
      qj <- c[c < q]
      qj1 <- c[q < c]
      if (length(qj) == 0) 
        qj <- 0
      h[x] <<- (q - qj[1])/(qj1[length(qj1)] - qj[1])
      qj5 <- inter$q5[c == qj[1]]
      qj51 <- inter$q5[c == qj1[length(qj1)]]
      (1 - h[x]) * qj5 + (h[x] * qj51)
    })
    
    q1 <- sapply(1:7, function(x) {
      q <- qx[x]
      c <- inter[[age_index[x]]]
      qj <- c[c < q]
      qj1 <- c[q < c]
      if (length(qj) == 0) 
        qj <- 0
      h[x] <<- (q - qj[1])/(qj1[length(qj1)] - qj[1])
      qj5 <- inter$q1[c == qj[1]]
      qj51 <- inter$q1[c == qj1[length(qj1)]]
      (1 - h[x]) * qj5 + (h[x] * qj51)
    })
    
  }
  
  data.frame(agegrp, women, child_born, child_dead, pi, di, 
             ki, qx, ti, year, q1, q5
             )
}




##  Exponencial: y(t) = y1 * exp( r * ( t2 - t1 ) ), 
##  r = log(y2/y1)/(t2-t1) es la tasa de crecimiento 
exp_interp <- 
  function( t, t1, t2, y1, y2 ){
    dt = t2 - t1
    r = log( y2 / y1 ) / dt
    y = y1 * exp( r * ( t - t1 ) )
    return( y )
  }