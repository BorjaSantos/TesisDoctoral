################################################################################################################################
#################################FUNCION 1: SIMULACION DE ESTUDIOS CON DATOS INDIVIDUALES#######################################
################################################################################################################################

# PARAMTEROS IMPORTANTES EN LA SIMULACION:
# 1. Numero de estudios
# 2. Numero de variables
# 3. Tamaño de la muestra
# 4. Media grupo tratamiento
# 5. SD grupo tratamiento
# 6. Media grupo control
# 7. SD grupo control.
# 8. Correlacion entre outcomes.

# VALORES DE LOS PARAMETROS
# Numero de estudios: 5, 10, 15, 25, 50.
# Numero de variables: 2, 4, 6.
# Tamano de la muestra: 20, 40, 60, 80, 100 (la mitad en cada rama) 
# Media grupo tratamiento distribuciones normales N(media,sd^2)
#   Var 1 (PANSS Total): N(60.56,16.97^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 30; Punt max = 210.
#   Var 2 (PANSS Positiva): N(13.65,5.90^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 7; Punt max= 49.
#   Var 3 (PANSS Negativa): N(17.82,6.92^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 7; Punt max= 49.
#   Var 4 (PANSS General): N(29.01,8.15^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 16; Punt max= 112.
#   Var 5 (HAMD): N(17.3,5.5) # Psychiatry Research 144 (2006) 57-63 Punt min = 0;Punt max= 52.
#   Var 6 (Calgary): N(17,15.3^2) # Banco de instrumentos para la psiquiatria clinica Punt min = 0; Punt max= 27.
# Media grupo control distribuciones normales N(media,sd^2)
# Parametros basados en las muestras de las publicaciones correspondientes necesarios para que no haya
# solapamiento entre los intervalos de confianza al 95% de los dos grupos.
#   Var 1 (PANSS Total): N(50.5,20.3^2) Punt min = 30; Punt max = 210.
#   Var 2 (PANSS Positiva): N(9,9.3^2) Punt min = 7; Punt max= 49.
#   Var 3 (PANSS Negativa): N(11.5,7^2) Punt min = 7; Punt max= 49.
#   Var 4 (PANSS General): N(20,9.15^2) Punt min = 16; Punt max= 112.
#   Var 5 (HAMD): N(8.4,4.6^2) Punt min = 0;Punt max= 52.
#   Var 6 (Calgary): N(3,4.2^2) Punt min = 0; Punt max= 27.
# Correlacion entre outcomes: 0, 0.1, 0.3, 0.5, 0.7, 0.9. 

# Libreri???as necesarias para la simulacion
library(truncdist) # De esta manera se truncaran los valores de las diferentes puntuaciones, para que no se simulen puntuaciones
# fuera de rango.
# Semilla para poder replicar los datos
set.seed(18052013)

# Funcion que va a simular los datos en el escenario que queramos

simulacion_datos <- function(n.estudios,n.vars,tamano.muestra,semilla,replicaciones){
  # Control de los paramtero de la funcion
  if(n.estudios > 50 || n.estudios < 5){
    stop("Numero de estudios incorrecto")
  }
  if((n.vars > 12 || n.vars < 4) && is.integer(n.vars/2)!= TRUE){
    stop("Numero de variables por estudio incorrecto")
  }
  if(tamano.muestra > 100 || tamano.muestra < 20){
    stop("Tamano de muestra incorrecto")
  }
  if(replicaciones < 5 || replicaciones > 150){
    stop("El numero de replicaciones es incorrecto")
  }
  # Funcion propiamente dicha
  database <- vector("list",replicaciones) # Lista donde se almacenaran las replicaciones del escenario deseado
  for(i in 1:replicaciones){
    database[[i]] <- vector("list",n.estudios) # Lista donde se almacenaran el número de estudios determinado
    set.seed(semilla)
    for(j in 1:n.estudios){
      Var1.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=30,b=210,mean=60.56,sd=16.97),0) # Puntuaciones PANSS total en tratamiento
      Var1.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=30,b=210,mean=50.5,sd=20.3),0) # Puntuaciones PANSS total en control
      Var2.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=13.65,sd=5.9),0) # Puntuaciones PANSS positiva tratamiento
      Var2.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=9,sd=9.3),0) # Puntuaciones PANSS positiva control
      Var3.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=17.82,sd=6.92),0) # Puntuaciones PANSS negativa tratamiento
      Var3.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=11,sd=5.7),0) # Puntuaciones PANSS negativa control
      Var4.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=16,b=112,mean=29.01,sd=8.15),0) # Puntuaciones PANSS general tratamiento
      Var4.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=16,b=112,mean=20.9,sd=9.15),0) # Puntuaciones PANSS general control
      Var5.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=52,mean=17.3,sd=5.5),0) # Puntuaciones HAMD tratamiento
      Var5.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=52,mean=8.4,sd=4.6),0) # Puntuaciones HAMD control
      Var6.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=27,mean=17.15,sd=15.3),0) # Puntuaciones Calgary tratamiento
      Var6.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=27,mean=3,sd=4.2),0) # Puntuaciones Calgary control
      database[[i]][[j]] <- as.data.frame(cbind(Var1.Trat,Var1.Ctrl,Var2.Trat,Var2.Ctrl,Var3.Trat,Var3.Ctrl,Var4.Trat,
                                                Var4.Ctrl,Var5.Trat,Var5.Ctrl,Var6.Trat,Var6.Ctrl))
      semilla <- semilla + 1
    }
  }
  # Para dimensionar en funcion del numero de variables requeridas
  for(i in 1:replicaciones){
    for(j in 1:n.estudios){
      if(n.vars == 4 || n.vars == 6 || n.vars == 8 || n.vars == 10 || n.vars == 12){
        database[[i]][[j]] <- database[[i]][[j]][,c(1:n.vars)]
      } else {
        stop("Numero de variables incorrecto")
      }
    }
  }
  return(database)
}

###############################################################################################################################
###################FUNCION 2: SIMULACION DE ESTUDIOS CON DATOS INDIVIDUALES (SIMPLIFICADA)#####################################
###############################################################################################################################

# COMPARACION DE LOS MÉTODOS DE METAANALISIS DESARROLLADOS FRENTE A LOS PROPUESTOS EN LA TESIS.
# LA COMPARACION SE HARA MEDIANTE SIMULACION BAJO DIFERENTES CONDICIONES.

# FECHA INICIO: 4 / MARZO / 2014
# FECHA FIN: 10/ MARZO / 2014

# PARAMTEROS IMPORTANTES EN LA SIMULACION:
# 1. Numero de estudios
# 2. Numero de variables
# 3. Tamaño de la muestra
# 4. Media grupo tratamiento
# 5. SD grupo tratamiento
# 6. Media grupo control
# 7. SD grupo control.
# 8. Correlacion entre outcomes.

# VALORES DE LOS PARAMETROS
# Numero de estudios: 5, 10, 15, 25, 50.
# Numero de variables: 2, 4, 6.
# Tamano de la muestra: 20, 40, 60, 80, 100 (la mitad en cada rama) 
# Media grupo tratamiento distribuciones normales N(media,sd^2)
#   Var 1 (PANSS Total): N(60.56,16.97^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 30; Punt max = 210.
#   Var 2 (PANSS Positiva): N(13.65,5.90^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 7; Punt max= 49.
#   Var 3 (PANSS Negativa): N(17.82,6.92^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 7; Punt max= 49.
#   Var 4 (PANSS General): N(29.01,8.15^2) # Rev Psiquiatr Salud Ment (Barc.) 2009;2(4):160-168 Punt min = 16; Punt max= 112.
#   Var 5 (HAMD): N(17.3,5.5) # Psychiatry Research 144 (2006) 57-63 Punt min = 0;Punt max= 52.
#   Var 6 (Calgary): N(17,15.3^2) # Banco de instrumentos para la psiquiatria clinica Punt min = 0; Punt max= 27.
# Media grupo control distribuciones normales N(media,sd^2)
# Parametros basados en las muestras de las publicaciones correspondientes necesarios para que no haya
# solapamiento entre los intervalos de confianza al 95% de los dos grupos.
#   Var 1 (PANSS Total): N(50.5,20.3^2) Punt min = 30; Punt max = 210.
#   Var 2 (PANSS Positiva): N(9,9.3^2) Punt min = 7; Punt max= 49.
#   Var 3 (PANSS Negativa): N(11.5,7^2) Punt min = 7; Punt max= 49.
#   Var 4 (PANSS General): N(20,9.15^2) Punt min = 16; Punt max= 112.
#   Var 5 (HAMD): N(8.4,4.6^2) Punt min = 0;Punt max= 52.
#   Var 6 (Calgary): N(3,4.2^2) Punt min = 0; Punt max= 27.
# Correlacion entre outcomes: 0, 0.1, 0.3, 0.5, 0.7, 0.9. 

# Librerias necesarias para la simulacion
library(truncdist) # De esta manera se truncaran los valores de las diferentes puntuaciones, para que no se simulen puntuaciones
# fuera de rango.
# Semilla para poder replicar los datos
set.seed(18052013)

# Función que va a simular los datos en el escenario que queramos

simulacion_datos2 <- function(n.estudios,n.vars,tamano.muestra,semilla,replicaciones){
  # Control de los paramtero de la funcion
  if(n.estudios > 50 || n.estudios < 5){
    stop("Numero de estudios incorrecto")
  }
  if((n.vars > 12 || n.vars < 4) && is.integer(n.vars/2)!= TRUE){
    stop("Numero de variables por estudio incorrecto")
  }
  if(tamano.muestra > 100 || tamano.muestra < 20){
    stop("Tamano de muestra incorrecto")
  }
  if(replicaciones < 5 || replicaciones > 150){
    stop("El numero de replicaciones es incorrecto")
  }
  # Funcion propiamente dicha
  database <- vector("list",replicaciones) # Lista donde se almacenaran las replicaciones del escenario deseado
  for(i in 1:replicaciones){
    database[[i]] <- vector("list",n.estudios) # Lista donde se almacenaran el número de estudios determinado
    set.seed(semilla)
    Var1.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=30,b=210,mean=60.56,sd=16.97),0) # Puntuaciones PANSS total en tratamiento
    Var1.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=30,b=210,mean=50.5,sd=20.3),0) # Puntuaciones PANSS total en control
    Var2.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=13.65,sd=5.9),0) # Puntuaciones PANSS positiva tratamiento
    Var2.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=9,sd=9.3),0) # Puntuaciones PANSS positiva control
    Var3.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=17.82,sd=6.92),0) # Puntuaciones PANSS negativa tratamiento
    Var3.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=7,b=49,mean=11,sd=5.7),0) # Puntuaciones PANSS negativa control
    Var4.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=16,b=112,mean=29.01,sd=8.15),0) # Puntuaciones PANSS general tratamiento
    Var4.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=16,b=112,mean=20.9,sd=9.15),0) # Puntuaciones PANSS general control
    Var5.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=52,mean=17.3,sd=5.5),0) # Puntuaciones HAMD tratamiento
    Var5.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=52,mean=8.4,sd=4.6),0) # Puntuaciones HAMD control
    Var6.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=27,mean=17.15,sd=15.3),0) # Puntuaciones Calgary tratamiento
    Var6.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=0,b=27,mean=3,sd=4.2),0) # Puntuaciones Calgary control
    database[[i]] <- as.data.frame(cbind(Var1.Trat,Var1.Ctrl,Var2.Trat,Var2.Ctrl,Var3.Trat,Var3.Ctrl,Var4.Trat,
                                         Var4.Ctrl,Var5.Trat,Var5.Ctrl,Var6.Trat,Var6.Ctrl))
    semilla <- semilla + (n.estudios) 
  }
  for(i in 1:replicaciones){
    if(n.vars == 4 || n.vars == 6 || n.vars == 8 || n.vars == 10 || n.vars == 12){
      database[[i]] <- database[[i]][,c(1:n.vars)]
    } else {
      stop("Numero de variables incorrecto")
    }
  }
  return(database)
}

##############################################################################################################################
####################FUNCION 3: g DE HEDGES' DE UNA SIMULACION CON REPLICACIONES (SIMPLIFICADA)################################
##############################################################################################################################

g_hedges_simulacion2_4vars <- function(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
                                       correlacion13=0,correlacion14=0,correlacion23=0,correlacion24=0,correlacion34=0){
  library(compute.es)
  n.vars <- 8
  datos <- simulacion_datos2(n.estudios,n.vars,tamano.muestra,semilla,replicaciones) 
  datos_meta <- diag(0,replicaciones,n.vars+6)
  for(i in 1:replicaciones){
    g1 <- mes(mean(datos[[i]][,1]),sd(datos[[i]][,1]),tamano.muestra,mean(datos[[i]][,2]),sd(datos[[i]][,2]),tamano.muestra)[12][,1]
    g2 <- mes(mean(datos[[i]][,3]),sd(datos[[i]][,3]),tamano.muestra,mean(datos[[i]][,4]),sd(datos[[i]][,4]),tamano.muestra)[12][,1]
    g3 <- mes(mean(datos[[i]][,5]),sd(datos[[i]][,5]),tamano.muestra,mean(datos[[i]][,6]),sd(datos[[i]][,6]),tamano.muestra)[12][,1]
    g4 <- mes(mean(datos[[i]][,7]),sd(datos[[i]][,7]),tamano.muestra,mean(datos[[i]][,8]),sd(datos[[i]][,8]),tamano.muestra)[12][,1]
    var.g1 <- mes(mean(datos[[i]][,1]),sd(datos[[i]][,1]),tamano.muestra,mean(datos[[i]][,2]),sd(datos[[i]][,2]),tamano.muestra)[13][,1]
    var.g2 <- mes(mean(datos[[i]][,3]),sd(datos[[i]][,3]),tamano.muestra,mean(datos[[i]][,4]),sd(datos[[i]][,4]),tamano.muestra)[13][,1]
    var.g3 <- mes(mean(datos[[i]][,5]),sd(datos[[i]][,5]),tamano.muestra,mean(datos[[i]][,6]),sd(datos[[i]][,6]),tamano.muestra)[13][,1]
    var.g4 <- mes(mean(datos[[i]][,7]),sd(datos[[i]][,7]),tamano.muestra,mean(datos[[i]][,8]),sd(datos[[i]][,8]),tamano.muestra)[13][,1]
    covar.g1g2 <- correlacion12*sqrt(var.g1)*sqrt(var.g2)
    covar.g1g3 <- correlacion13*sqrt(var.g1)*sqrt(var.g3)
    covar.g1g4 <- correlacion14*sqrt(var.g1)*sqrt(var.g4)
    covar.g2g3 <- correlacion23*sqrt(var.g2)*sqrt(var.g3)
    covar.g2g4 <- correlacion24*sqrt(var.g2)*sqrt(var.g4)
    covar.g3g4 <- correlacion34*sqrt(var.g3)*sqrt(var.g4)
    input <- c(g1,g2,g3,g4,var.g1,covar.g1g2,covar.g1g3,covar.g1g4,var.g2,covar.g2g3,covar.g2g4,var.g3,covar.g3g4,var.g4)
    datos_meta[i,] <- input
  }
  colnames(datos_meta) <- c("g1","g2","g3","g4","var.g1","covar.g1g2","covar.g1g3","covar.g1g4","var.g2","covar.g2g3",
                            "covar.g2g4","var.g3","var.g3g4","var.g4")
  return(as.data.frame(datos_meta))
}

##############################################################################################################################
##########################FUNCION 4: g DE HEDGES' DE UNA SIMULACION CON REPLICACIONES#########################################
##############################################################################################################################

g_hedges_simulacion_4vars <- function(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
                                      correlacion13=0,correlacion14=0,correlacion23=0,correlacion24=0,correlacion34=0){
  library(compute.es)
  n.vars <- 8
  datos <- simulacion_datos(n.estudios,n.vars,tamano.muestra,semilla,replicaciones)
  datos_meta <- vector("list",replicaciones)
  g1 <- vector(mode="numeric",length=n.estudios)
  g2 <- vector(mode="numeric",length=n.estudios)
  g3 <- vector(mode="numeric",length=n.estudios)
  g4 <- vector(mode="numeric",length=n.estudios)
  var.g1 <- vector(mode="numeric",length=n.estudios)
  var.g2 <- vector(mode="numeric",length=n.estudios)
  var.g3 <- vector(mode="numeric",length=n.estudios)
  var.g4 <- vector(mode="numeric",length=n.estudios)
  covar.g1g2 <- vector(mode="numeric",length=n.estudios) 
  covar.g1g3 <- vector(mode="numeric",length=n.estudios)
  covar.g1g4 <- vector(mode="numeric",length=n.estudios) 
  covar.g2g3 <- vector(mode="numeric",length=n.estudios)
  covar.g2g4 <- vector(mode="numeric",length=n.estudios) 
  covar.g3g4 <- vector(mode="numeric",length=n.estudios) 
  for(i in 1:replicaciones){
    for(j in 1:n.estudios){
      g1[j] <- mes(mean(datos[[i]][[j]][,1]),sd(datos[[i]][[j]][,1]),tamano.muestra,mean(datos[[i]][[j]][,2]),sd(datos[[i]][[j]][,2]),tamano.muestra)[12][,1]
      g2[j] <- mes(mean(datos[[i]][[j]][,3]),sd(datos[[i]][[j]][,3]),tamano.muestra,mean(datos[[i]][[j]][,4]),sd(datos[[i]][[j]][,4]),tamano.muestra)[12][,1]
      g3[j] <- mes(mean(datos[[i]][[j]][,5]),sd(datos[[i]][[j]][,5]),tamano.muestra,mean(datos[[i]][[j]][,6]),sd(datos[[i]][[j]][,6]),tamano.muestra)[12][,1]
      g4[j] <- mes(mean(datos[[i]][[j]][,7]),sd(datos[[i]][[j]][,7]),tamano.muestra,mean(datos[[i]][[j]][,8]),sd(datos[[i]][[j]][,8]),tamano.muestra)[12][,1]
      var.g1[j] <- mes(mean(datos[[i]][[j]][,1]),sd(datos[[i]][[j]][,1]),tamano.muestra,mean(datos[[i]][[j]][,2]),sd(datos[[i]][[j]][,2]),tamano.muestra)[13][,1]
      var.g2[j] <- mes(mean(datos[[i]][[j]][,3]),sd(datos[[i]][[j]][,3]),tamano.muestra,mean(datos[[i]][[j]][,4]),sd(datos[[i]][[j]][,4]),tamano.muestra)[13][,1]
      var.g3[j] <- mes(mean(datos[[i]][[j]][,5]),sd(datos[[i]][[j]][,5]),tamano.muestra,mean(datos[[i]][[j]][,6]),sd(datos[[i]][[j]][,6]),tamano.muestra)[13][,1]
      var.g4[j] <- mes(mean(datos[[i]][[j]][,7]),sd(datos[[i]][[j]][,7]),tamano.muestra,mean(datos[[i]][[j]][,8]),sd(datos[[i]][[j]][,8]),tamano.muestra)[13][,1]
      covar.g1g2[j] <- correlacion12*sqrt(var.g1[j])*sqrt(var.g2[j])
      covar.g1g3[j] <- correlacion13*sqrt(var.g1[j])*sqrt(var.g3[j])
      covar.g1g4[j] <- correlacion14*sqrt(var.g1[j])*sqrt(var.g4[j])
      covar.g2g3[j] <- correlacion23*sqrt(var.g2[j])*sqrt(var.g3[j])
      covar.g2g4[j] <- correlacion24*sqrt(var.g2[j])*sqrt(var.g4[j])
      covar.g3g4[j] <- correlacion34*sqrt(var.g3[j])*sqrt(var.g4[j])
      input <- cbind(g1,g2,g3,g4,var.g1,covar.g1g2,covar.g1g3,covar.g1g4,var.g2,covar.g2g3,covar.g2g4,var.g3,covar.g3g4,var.g4)
    }
    datos_meta[[i]] <- as.data.frame(input)
  }
  return(datos_meta)
}

##############################################################################################################################
##########################FUNCION 5: METAANALISIS MULTIVARIADO DE SIMULACION (MODIFICADA)#####################################
##############################################################################################################################

###############
# FUNCION 5.1 #
###############

#datos <- g_hedges_simulacion2_4vars(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
#                                    correlacion13=0,correlacion14=0,correlacion23=0,correlacion24=0,correlacion34=0)

meta_multi_simulacion2_4vars <- function(data,metodo="reml"){
  library(mvmeta)
  meta_multi <- mvmeta(cbind(g1,g2,g3,g4),as.data.frame(data)[5:14],data=data,method=metodo)
  coeficientes <- meta_multi$coefficients
  coeficientes_inferencia <- summary(meta_multi)$coefficients
  coeficientes_var_cov <- summary(meta_multi)$corRandom
  output <- list(coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  return(output)
}

###############
# FUNCION 5.2 #
###############

meta_multi_simulacion2_4vars_bis <- function(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
                                             correlacion13=0,correlacion14=0,correlacion23=0,correlacion24=0,correlacion34=0,
                                             metodo="reml"){
  library(mvmeta)
  datos <- g_hedges_simulacion2_4vars(n.estudios,tamano.muestra,semilla,replicaciones,correlacion12,
                                      correlacion13,correlacion14,correlacion23,correlacion24,correlacion34)
  data <- as.data.frame(datos)
  meta_multi <- mvmeta(cbind(g1,g2,g3,g4),as.data.frame(data)[5:14],data=data,method=metodo)
  coeficientes <- meta_multi$coefficients
  coeficientes_inferencia <- summary(meta_multi)$coefficients
  coeficientes_var_cov <- summary(meta_multi)$corRandom
  output <- list(coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  return(output)  
}

##############################################################################################################################
##########################FUNCION 6: METAANALISIS MULTIVARIADO DE SIMULACION (UNO POR REPLICACION)############################
##############################################################################################################################

###############
# FUNCION 6.1 #
###############

#datos <- g_hedges_simulacion_4vars(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
#                                   correlacion13=0,correlacion14=0,correlacion23=0,correlacion24=0,correlacion34=0)

meta_multi_simulacion_4vars <- function(data,metodo="reml",replicaciones=5){
  library(mvmeta)
  meta_multi <- vector(mode="list",length=replicaciones)
  for(i in 1:length(data)){
    datos <- data[[i]]
    meta_resultados <- mvmeta(cbind(g1,g2,g3,g4),S=as.data.frame(data)[5:14],data=datos,method=metodo)
    coeficientes <- meta_resultados$coefficients
    coeficientes_inferencia <- summary(meta_resultados)$coefficients
    coeficientes_var_cov <- summary(meta_resultados)$corRandom
    meta_multi[[i]] <- list(dat=data[[i]],coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  }
  return(meta_multi)  
}

###############
# FUNCION 6.2 #
###############

meta_multi_simulacion_4_vars_bis <- function(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
                                             correlacion13=0,correlacion14=0,correlacion23=0,correlacion24=0,correlacion34=0,
                                             metodo="reml"){
  library(mvmeta)
  datos <- g_hedges_simulacion_4vars(n.estudios,tamano.muestra,semilla,replicaciones,correlacion12,correlacion13,correlacion14,
                                     correlacion23,correlacion24,correlacion34)
  meta_multi <- vector(mode="list",length=replicaciones)
  for(i in 1:replicaciones){
    data <- as.data.frame(datos[[i]])
    meta_resultados <- mvmeta(cbind(g1,g2,g3,g4),S=as.data.frame(data)[5:14],data=data,method=metodo)
    coeficientes <- meta_resultados$coefficients
    coeficientes_inferencia <- summary(meta_resultados)$coefficients
    coeficientes_var_cov <- summary(meta_resultados)$corRandom
    meta_multi[[i]] <- list(dat=datos[[i]],coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  }
  return(meta_multi)  
}