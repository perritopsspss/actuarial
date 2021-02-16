TablaSBS = read.table("TMSPP2017.txt",T)
attach(TablaSBS)
library(lifecontingencies)

#mujeres sanas para el 2020
probas = as.vector(Tabla$SPPS2017M*(1-Tabla$AaxM)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla2 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla2")


#hombre invalido
probas = as.vector(Tabla$SPPI2017H*(1-Tabla$AaxH)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla3 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla3")


#mujer invalida
probas = as.vector(Tabla$SPPI2017M*(1-Tabla$AaxM)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla4 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla4")

#hombre sano 
probas = as.vector(Tabla$SPPH2017H*(1-Tabla$AaxH)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla3 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla3")

###############################################################################################

