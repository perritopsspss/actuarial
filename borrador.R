#Considere el caso de 2 hermanos varones saludables de 18 y 25 años en el 2021
#Si contratan un seguro de tal modo que al fallecer uno de ellos, el otro cobra 
#una indemnización de 200 mil soles y el rendimiento anual es del 2.5%, la prima
#única a ser pagada es: 

TablaSBS = read.table("TMSPP2017.txt",T)
attach(TablaSBS)
library(lifecontingencies)
Tabla1 = probs2lifetable(probs = SPPS2017H*(1-AaxH)^4,
                         radix = 10^6, type = "qx")
Tabla2 = probs2lifetable(probs = SPPS2017H*(1-AaxH)^4,
                         radix = 10^6, type = "qx")

i=2.5/100
Tablas3 = list(T5 = Tabla1, T6 = Tabla2)
A = Axyzn(Tablas3, x = c(18,25), i = i, status = "joint")
A  #esta es la prima unitaria 

S = 200000
P = S*A
P #monto de la prima única



i = 1.25/100
d = i/(1+i)
A = 0.032+0.875
a = (1-0.940)/d
A/a

40000/A



#En el año 2021, una mujer sana de 34 años contrata un seguro para que sus beneficiarios reciban____
#soles si fallece teniendo al menos 64 pero menos de 80 años. La TEA considerada es del 2% 
#y se realizará el pago de 10 primas anuales (como máximo) de 1115.87 soles cada una.
#Use las tablas de la SBS

> Tabla  = read.table("TMSPP2017.txt",T)
> probas = as.vector(Tabla$SPPS2017M*(1-Tabla$AaxM)^4)
> Tabla1 = probs2lifetable(probs = probas, 
                           +                          radix = 10^6,  
                           +                          type  = "qx",  
                           +                          name = "Tabla1")
> (A = Axn(Tabla1, x = 34, m = 30, n = 16, i = 0.02))
[1] 0.06800987
> (a = axn(Tabla1, x = 34, n = 10, i = 0.02))
[1] 9.142183
> A/a # Prima unitaria anual
[1] 0.007439128
> 1115.87*a/A
[1] 150000.1



#################################################################


#mujer invalida
probas = as.vector(TablaSBS$SPPI2017M*(1-TablaSBS$AaxM)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla4 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                        name = "Tabla4")
#17q31
qxt(Tabla4, 31, 17)

#m|nQx = mPx * nQx+m

#32|1q31 
x= 31
m=32
n=1

pxt(Tabla4, x = x, t = m) * qxt(Tabla4, x = x+m, t = n)
pxt(Tabla4,31,32)* qxt(Tabla4,63 ,1)
###################################################################


#################################################################
probas = as.vector(TablaSBS$SPPS2017H*(1-TablaSBS$AaxH)^4)
Tabla3 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla3")

i = 3/100
d=  i/(1+i)
(P = axn(Tabla3, x = 25, i = i))

(P = (1-Axn(Tabla3, x=25, i=i))/d)

18000*P

#########################################################
#10p30 =exp(- integral de 0 hasta 10 0.6ds) = exp(-0.6(10-0)) = exp(-6) = 0.0025

mu = function(s) return(0.3+0*s)
(I = integrate(mu,0,5)$value)
exp(-I)

########################################################
#esposa sana 
probas = as.vector(TablaSBS$SPPS2017M*(1-TablaSBS$AaxM)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla2 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla2")

#hombre sano 

probas = as.vector(TablaSBS$SPPS2017H*(1-TablaSBS$AaxH)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla3 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla3")

Tablas = list(T1 = Tabla2, T2 = Tabla3)

x = c(45, 44)
t = 14
npx  = pxt(Tablas[[1]], x=x[1], t=t)
npy  = pxt(Tablas[[2]], x=x[2], t=t)
npxy = pxyzt(Tablas, x=x, t=t, status = "joint")
npx + npy - npxy

pxyzt(Tablas, x=x, t=t, status = "last")

i = 0.03
A = Axyzn(Tablas, x = x, i = i, status = "joint")
A  #esta es la prima unitaria 

S = 300000
P = S*A
P #monto de la prima única

######################################################
100000*0.85



#############################3



probas17 = as.vector(TablaSBS$SPPS2017M*(1-TablaSBS$AaxM)^0)
probas18 = as.vector(TablaSBS$SPPS2017M*(1-TablaSBS$AaxM)^1)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla2_17 = probs2lifetable(probs = probas17,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla2")

Tabla2_18 = probs2lifetable(probs = probas18,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla2")


q17 = qxt(Tabla2_17, 25, 1)
q18 = qxt(Tabla2_18, 25, 1)

################################################
#pregunta 8 

edad = 10.12
Sx = function(x) (1 - (x^0.25)/(105^0.25) )^0.1
f = function(m) return(Sx(edad+m)/Sx(edad)-0.5)
mediana = uniroot(f, lower = 0, upper = 105-edad)$root


num = integrate(Sx, lower = edad, upper = 105)$value
den = Sx(edad)
(e30 = num/den)
