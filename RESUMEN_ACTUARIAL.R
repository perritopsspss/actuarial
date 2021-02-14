
# ACTUARIAL_TEMAS_ EXAMEN FINAL_RECUPERACION ------------------------------


# perrito habla huevadas

# TEMA: ANÁLISIS DE TIEMPOS DE VIDA ---------------------------------------
Tabla1 = read.table("TMSPP2017.txt",T)
library(lifecontingencies)
probas = as.vector(Tabla1$SPPS2017H*(1-Tabla1$AaxH)^3)
Tabla1 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla1")

dx = rep(0,length(Tabla@lx))
for(i in 1:(length(Tabla@lx))){dx[i] = dxt(Tabla1,i-1,1)}
library(psych)
headTail(cbind(Tabla1@lx,dx))

(d0 = dx[1])
dxt(Tabla1,0,1)
# 7523.007
#d0 = De un total de 1 millón de varones sanos, se espera que 7523 de ellos
#fallezcan antes de cumplir un año en el 2020.



# Px y Qx -----------------------------------------------------------------

#Px = proporcion de individuos de la edad x que alcanza la edad x+1}
#Qx es la proporción de defunciones de los individuos de edad x que no alcanzan la edad x + 1:

#pxt(Tabla1,0,1) =  lx[2]/lx[1]
#proba de que en el 2020, un varon sano llegue a vivir al menos un año más 

#pxt(Tabla1,45,1) = lx[47]/lx[46] = 0.9984188
#proba de que en el 2020, un varon sano de 45 años llegue a vivir al menos 
#un año más o que fallezca despues de haber cumplido 46 años 

#qxt(Tabla1,0,1) =  1-lx[2]/lx[1] = 0.007523007
#q0 = 0, proba de que en el 2020, un varon sano fallezca antes de cumplir un año más



# nPx y nQx ---------------------------------------------------------------

#nPx es la probabilidad de que una persona de edad x viva al menos n años más

pxt(Tabla1, 43, 30)


#nQx es la probabilidad de que una persona de edad x fallezca antes de alcanzar la
#edad x + n.

#               nQx = 1 − nPx

qxt(Tabla1, 28, 37)


# m|nQx  ----------------------------------------------------------------

#m|nQx, probabilidad de que un sujeto de edad x alcance la edad de x + m años
#y fallezca antes de los n siguientes (antes de la edad x + m + n).

#               m|nQx = mPx * nQx+m

#  m−1|qx, proba de que un sujeto de edad x fallezca antes de alcanzar la edad x+m



# K(x) --------------------------------------------------------------------

#es el tiempo de vida futuro truncado o numero de años completos antes de morir

#                         K(x) = [X-x], siendo x la edad actual del asegurado

#si el asegurado tiene x = 41 años y fallece a los X = 44.25 años,
#entonces k = K(41) = [44.25 -41] = 3 

#                        P(K(x) = k) = k|qx

#ejemplo:

#P(K(60) = 24) = 24|q60
pxt(Tabla1,60,24) * qxt(Tabla1,84,1)

#P(K(25) <= 45) = 46q25
qxt(Tabla1, 25,46)

#SK(50)(30) = P(K(50) > 30) = 31p50
pxt(Tabla1,50,31)

#P(12 <=  K(80) < 25) = 12|13q80
pxt(Tabla1,80,12)*qxt(Tabla1,92,13)







# ex ----------------------------------------------------------------------

#la esperanza de K(x) es conocida como la esperanza de vida restante 

#            ex = sumatoria de npx, del 1 al infinito 

#e50 = 34.59, es decir se espera que un varón sano de 50 años viva 34.59 años más.

exn(Tabla1,x=50)

#e0 = 81.65, es decir se espera de vida al nacer de un varón sano es 81.65 años

exn(Tabla1,x=0)




# FUNCIONES CONTINUAS DE SUPERVIVENCIA ------------------------------------

#                        S(x) = 1 − (x^2/11025)   para para 0 <= x <= 105

Sx = function(x){return(1-x^2/11025)}
x = seq(0,105,0.01)
S = Sx(x)
d = data.frame(x,S)



# Px y Qx -----------------------------------------------------------------

#px = S(x+1)/S(x)

#p45 = S(46)/S(45)
Sx(46)/Sx(45)

#qx = 1 - px = 1 - (S(x+1)/S(x))
1-Sx(64)/Sx(63)















# nPx y nQx ---------------------------------------------------------------

#        nPx = S(x+n)/S(x)

# 20p10
(p.2010 = Sx(30)/Sx(10))

#        nqx = 1 -  nPx

# 4q10 = 1 - 4p10
(q.410 = 1-Sx(14)/Sx(10))




# m|nQx -------------------------------------------------------------------

# m|nqx = S(x + m) − S(x + m + n)/(S(x))
# 1|4q0

x = 0
m= 1
n= 4
prob = (Sx(x+m)-Sx(x+m+n))/Sx(x)
prob


# T(x) --------------------------------------------------------------------

#es el tiempo de vida restante de la persona

#           T(x) = X - x


# ex ----------------------------------------------------------------------

Sx = function(x){1-x^2/11025}
num = integrate(Sx, lower = 30, upper = 105)$value
den = Sx(30)
(e30 = num/den)

#las personas de 30 años vivirán en promedio 44.44 años más. La esperanza de
#vida de las personas de 30 años es de 74.44 años.


# mx ----------------------------------------------------------------------

#el tiempo de vida futuro mediano a la edad de x años,
#se obtiene al resolver la ecuacion: 

#P(T(x) > mx ) = 0.5
#S(x + mx)/S(x) = 0.5

#Cálculo de la mediana del tiempo de vida restante a lo 50 años
#m50:

S(50 + m50)/S(50) = 0.5

Sy = function(y){1-y^2/11025}
f = function(m) return(Sy(50+m)/Sy(50)-0.5)
uniroot(f, lower = 0, upper = 105-50)$root

#corroborando:
Sy(82.23442)/Sy(50)
#0.5



# Ux ----------------------------------------------------------------------
#fuerza de mortalidad 

#Ux = F(x)/S(x)

#esta función guarda relacion con tpx y tqx

#   tpx = exp(- integral de 0 a t Ux + sds)

#Si Ux = 0.6, habllar el valor de 10p30

#10p30 =exp(- integral de 0 hasta 10 0.6ds) = exp(-0.6(10-0)) = exp(-6) = 0.0025

mu = function(s){return(0.6+0*s)}
(I = integrate(mu,0,10)$value)
exp(-I)

#Si Ux = X/100, hallar el valor de 7q42

#7q42 = exp( - integral de 0 al 7 (42 + s/100) = exp (-((42*7)/100 + 7^2/200)) = exp ( -3.185) = 0.0414

mu = function(s){return((42+s)/100)}
(I = integrate(mu,0,7)$value)

1-exp(-I)




# TEMA: MODELOS DE SEGUROS DE VIDA  ---------------------------------------------


# SEGURO_DOTAL_PURO -------------------------------------------------------

#(v^n)(npx )

#CASO DISCRETO

TablaSBS = read.table("TMSPP2017.txt",T)
attach(TablaSBS)
library(lifecontingencies)
Tabla1 = probs2lifetable(probs = SPPS2017H*(1-AaxH)^3,
                         radix = 10^6, type = "qx")
#x = edad
#n = los años que faltan para que la persona de x años reciba el seguro
(EZ = Exn(Tabla1, x=22, n=50, i=0.036))

#si quiere recibir 50 mil soles:
EZ*50000

#CASO CONTINUO  

#Determine el monto de la prima que debería pagar una persona de 40 años para recibir
#80 mil soles si se encuentra viva a los 70 años. Considere una tasa efectiva anual de
#4.05%

#S(x) = exp(−0.0003 × (exp(0.1x) − 1))

S = function(x) exp(-0.0003*(exp(0.1*x)-1))
v = 1/1.0405
p = S(70)/S(40)
(E = v^30*p)



# SEGURO_VITALICIO --------------------------------------------------------

#caso discreto 

#Determine la prima que debe pagar una señora inválida de 52 años para que cuando fallezca,
#su familia esté asegurada por un monto de 200 mil soles, considerando una tasa de 3% anual.

# x = 52 
Tabla2 = probs2lifetable(probs = SPPI2017M*(1-AaxM)^3, radix = 10^6, type = "qx")
(EZ = Axn(Tabla2, x=52, i=0.03))

EY = 200000*EZ

#caso continuo 
#E(Zx) = Ax = E (exp (−deltaT(x)))
#SX (x) = 1 − (x^3/10^6)

#Determine la prima que debe pagar una persona de 60 años para asegurar a su familia
#en caso de su muerte por un monto de 200 mil soles. Considere una TEA del 2.6%

integrando = function(t) exp(-log(1.026)*t)*3*(60+t)^2/(10^6-60^3)
x = 60 #edad
w = 100 # hasta donde llega el x en la función 
(A60 = integrate(integrando,0,w-x)$value)

200000*A60



# SEGURO_TEMPORAL ---------------------------------------------------------

#Caso discreto 

#Determine la prima que debe pagar una señora sana de 36 años para asegurar a su familia por
#un monto de 95 mil soles siempre y cuando fallezca antes de cumplir 100 años. Considere una tasa de
#4% anual.

Tabla3 = probs2lifetable(probs = SPPS2017M*(1-AaxM)^3, radix = 10^6, type = "qx")
# x = edad, n= años que le falta a la persona para terminar el tiempo que dura su seguro
(EZ = Axn(Tabla3, x=36, n=64, i=0.04))
(EY = 95000*EZ)

#Caso continuo 
#E(Zx) = Ax = E (exp (−deltaT(x)))
#SX (x) = 1 − (x^3/10^6)

#Determine la prima que debe pagar una persona de 63 años para asegurar a su familia
#por un monto de 110 mil soles siempre y cuando fallezca dentro de los próximos 25 años.
#Considere una tasa de 3% anual.

integrando = function(t) exp(-log(1+0.03)*t)*3*(60+t)^2/(10^6-60^3)
x = 63 #edad
w = 88 # 25+63
(A60 = integrate(integrando,0,w-x)$value)
200000*A60


# SEGURO_DOTAL_MIXTO ------------------------------------------------------

#Caso discreto 

#seguro temporal + seguro dotal puro 
#Ax:n + nEx
#Determine la prima que debe pagar un señor inválido de 59 años por un seguro
#dotal (mixto) de 20 años y un monto de 99 mil soles. Considere una tasa de 3.5% anual.

Tabla4 = probs2lifetable(probs = SPPI2017H*(1-AaxH)^3,
                         radix = 10^6, type = "qx")
# x=edad , n=20 es el plazo que dura el temporal 
(EZ = AExn(Tabla4, x=59, n=20, i=0.035))
(EY = 99000*EZ)

#otra forma 

(TEMP = Axn(Tabla4, x=59, n=20, i=0.035))
#[1] 0.4162387
(DOTP = Exn(Tabla4, x=59, n=20, i=0.035))
#[1] 0.1954417
(DOTM = TEMP + DOTP)
#[1] 0.6116803

#Caso continuo 

#Determine la prima que debe pagar una persona de 58 años que contrata un seguro de
#vida dotal de 20 años por un monto de 100 mil soles, a una tasa de interés nominal
#anual de 3% capitalizable semestralmente.

#La tasa nominal anual de 3% capitalizable semestralmente = 1.5% semestral = 3.0225% anual.

#primero el temporal 
integrando = function(t){
  exp(-log(1.03225)*t)*1/(4*42)*(1-t/42)^-0.75}
x = 58 #edad
n = 20 #la duración del seguro dotal
(A58 = integrate(integrando,0,n)$value)

#segundo el dotal puro 
#  20E58 = 20P58/(1 + 0.03225)^20

Sx = function(x){
  return((1-x/100)^0.25)
} 

(E58 = (Sx(78)/Sx(58))/(1+0.03225)^20)

A58+E58

# SEGURO_VITALICIO_DIFERIDO ---------------------------------------------------------

#Seguro dotal puro +  seguro vitalicio

#Caso discreto

#Determine la prima que debe pagar una señora sana de 55 años por un monto
#asegurado de 80 mil soles, si éste será efectivo desde que cumpla 65 años. Considere
#una tasa de 3.2% anual.

#x = 55, m = 10 (años que le faltan para que se haga efectivo el seguro)
Tabla5 = probs2lifetable(probs = SPPS2017M*(1-AaxM)^3,
                         radix = 10^6, type = "qx")
(EZ = Axn(Tabla5, x=55, m=10, i=0.032))
(EY = 80000*EZ)

#de otra forma 
#m|Ax =mEx × Ax+m
#10|A55 =10E55 × A65

(E = Exn(Tabla5, x=55, n=10, i=0.032))
# 0.7057307
(A = Axn(Tabla5, x=65, i=0.032))
# 0.4719536
E*A

#Caso continuo 

SX (x) =  (1 − x/100)^1/4
#Determine la prima que debe pagar una persona de 62 años que
#contrata un seguro de vida completa diferido 16 años por un monto de 100 mil soles, a
#una tasa de interés nominal bimestral de 0.4% capitalizable mensualmente.

#La tasa nominal bimestral de 0.4% capitalizable mensualmente = 0.2% mensual = 2.4266% anual.

integrando = function(t){
  exp(-log(1.024266)*t)*1/(4*38)*(1-t/38)^-0.75}
x = 62
w = 100
(A60 = integrate(integrando,16,w-x)$value)
100000*A60

#de otra manera
#16|A62 =16E62 × A78

Sx = function(x){return((1-x/100)^0.25)}
(E62 = (Sx(78)/Sx(62))/(1+0.024266)^16)
## [1] 0.5943683
integrando = function(t){
  exp(-log(1.024266)*t)*1/(4*22)*(1-t/22)^-0.75}
(A78 = integrate(integrando,0,100-78)$value)
## [1] 0.6626973
E62*A78
## [1] 0.3938863



# SEGURO_DOTAL_DIFERIDO ---------------------------------------------------

#Caso discreto

#Determine la prima que debe pagar una señora inválida de 55 años por un
#monto asegurado de 80 mil soles, si éste será efectivo siempre y cuando fallezca mientras
#tenga entre 76 y 88 años. Considere una tasa de 2.5% anual.

Tabla2 = probs2lifetable(probs = SPPI2017M*(1-AaxM)^3, radix = 10^6, type = "qx")
(E = Exn(Tabla2, x=55, n= 76-55, i=0.032))
(EZ = Axn(Tabla2, x=55+21, n=88-(55+21), i=0.04))
E*EZ

#Caso continuo 
#(1-x/100)^0.25

Sx = function(x){return((1-x/100)^0.25)}
(E55 = (Sx(76)/Sx(55))/(1+0.025)^21)

integrando = function(t)
  {exp(-log(1.025)*t)*1/(4*24)*(1-t/24)^-0.75}
x = 55+21
n = 88-76
(A55 = integrate(integrando,0,n)$value)

A55*E55





# TEMA: ANUALIDADES



# TEMA: ANUALIDADES -------------------------------------------------------

# ANUALIDADES_VITALICIAS_ANTICIPADAS --------------------------------------------------
# La probabilidad de pago en cada posible fecha ---------------------------

#Sumatoria de 1/(1 + i)^t × tpx (desde T al w-x)= tEx
#para x = 40: ¨a40 = 1 + 1/(1 + i) × p40 + ... + 1/(1 + i)70 × 70p40

#x = 108: 1 + E108 + 2E108

#El valor presente de la anualidad anticipada cierta:
#v=1/(1+i)
#d=i/(1+i)
#¨an| = (1 − V^n)d

#Suponga que está sujeto a morir según la tabla de mortalidad de la SBS.
#ejercicio:
#Un hombre sano de 55 años desea contratar un seguro de vida completa a fin
#de recibir anualidades de 18000 soles al inicio de cada año, considerando una TEA de
#2%, calcule el valor de la prima bajo las siguientes condiciones:
#a) Suponga que vivirá 110 años.

#¨an| = (1 − V^n)d
i = 0.02
v = 1/(1+i)
d = i/(1+i)
n = 55
(a = (1 - v^n)/d)
#33.83828

#¨a55|(con una rayita encima de la n) = 33.83
#33.83 es el  monto que hay que pagar para recibir 55 pagos anticipados de 1 sol.

# La definición de valor presente cierto de una anualidad de pagos anticipados --------

# ¨a = 1/d (1 − Ax )
#Ax es el valor presente actuarial de un seguro de vida completa
library(lifecontingencies)

TablaSBS = read.table("TMSPP2017.txt",T)
attach(TablaSBS)
Tabla = probs2lifetable(probs = SPPS2017H*(1-AaxH)^4, radix = 10^6, type = "qx")
(P = axn(Tabla, x = 55, i = 0.02))

#otra forma:
i = 0.02
d = i/(1+i)
(P = (1-Axn(Tabla,x=55,i=0.02))/d)

1800*P

#interpretacion

#a55 (con dos puntos sobre la a) = 22.81:
#22.81 es el monto que hoy una persona de 55 años debe pagar para recibir pagos 
#anticipados vitalicios de un sol anual a partir de los 61 años.



# ANUALIDADES_VITALICIAS_VENCIDAS -----------------------------------------

#El valor presente de una anualidad de n pagos vencidos:
#v=1/(1+i)
#d=i/(1+i)
#an| = (1 − V^n)i

#El valor presente actuarial de la anualidad:
#¨ax-1

#ejercicio:
#Una mujer inválida de 60 años desea contratar un seguro de vida completa a
#fin de recibir pagos de 15000 soles al final de cada año, condicionales a su sobrevivencia.
#Considerando una TEA de 1.4%, calcule el valor de la prima bajo las siguientes
#condiciones:

# Considerando la probabilidad de muerte año a año ------------------------
i = 0.014
v = 1/(1+i)
d = i/(1+i)
n = 50
(a = (1 - v^n)/i)
#35.78557
a*15000

#interpretación:
#a50| (con rayita encima): 35.79 es el monto que hay que pagar para recibir 50 pagos vencidos de 1 sol. 

# Suponga que está sujeta a morir según la tabla de mortalidad de la SBS  --------

library(lifecontingencies)
TablaSBS = read.table("TMSPP2017.txt",T)
attach(TablaSBS)
Tabla = probs2lifetable(probs = SPPI2017M*(1-AaxM)^4, radix = 10^6, type = "qx")
(P = axn(Tabla, x=60, i=0.014) - 1)

#otra forma:
(P = axn(Tabla, x=60, i=0.014, payment = "arrears"))
#otra forma:
(P = (1-(1+i)*Axn(Tabla,x=60,i=0.014))/i)

#Interpretacion: 
# a60 = 16.81 es el monto que hay que pagar para recibir pagos vencidos vitalicios
#de un sol anual a partir de los 61 años
# ANUALIDADES_VITALICIAS_CONTINUAS ----------------------------------------
#El valor presente actuarial de la anualidad continua:

# ax (con raya encima de la a) = (1 − Ax)/gama
# gama = log (1 + i)

#ejercicio:

#Calcule el valor de la prima que debe pagar una persona de 43 años cuyo
#tiempo de vida puede ser explicado mediante la función de supervivencia

#S(x) = 1 − (x^2.1)/17209.57

#para 0 <= x <= 104, a fin de recibir 20 mil soles anuales (de manera continua, por ejemplo
#500 soles semanales) mientras se encuentre vivo. Considere una TEA de 2%.

S = function(x){1-(x/104)^2.1}
d = log(1+0.02)
f = function(t){exp(-d*t)*S(43+t)/S(43)}
(P = integrate(f,0,104-43)$value)
#23.80
P*20000

#IMPORTANTE:
# a43 (con raya encima de la a) = (1 − A43)/d




# ANUALIDADES_TEMPORALES_ANTICIPADAS --------------------------------------------------

# ¨ax:n|

# Considerando la definición de valor presente cierto de una anualidad de pagos anticipados --------

#   (1 − Ax:n|)/d

#Un hombre sano de 52 años desea contratar un seguro de vida temporal de 30
#años a fin de recibir anualidades de 15000 soles al inicio de cada año, considerando una
#TEA de 3%, calcule el valor de la prima bajo las siguientes condiciones


# Suponga que el hombre vivirá, con certeza, 82 años. ---------------------

i = 0.03
v = 1/(1+i)
d = i/(1+i)
n = 30 #tiempo que debe transcurrir para que se haga efectivo sus anualidades (es temporal)
(a = (1 - v^n)/d)
a*15000

#interpretacion
# ¨a30|(falta la rayita encima del 30) = 20.19
# 20.19 es el  monto que hay que pagar para recibir 30 pagos temporales anticipados de 1 sol.

# Suponga que el hombre está sujeto a morir según la tabla de mortalidad --------
Tabla = probs2lifetable(probs = SPPS2017H*(1-AaxH)^4, radix = 10^6, type = "qx")
(P = axn(Tabla, x=52, n=30, i=0.03))

#otra forma:
(P = (1-AExn(Tabla, x=52, n=30, i=0.03))/d)

#Interpretacion
# ¨a52:30| (con rayita encima del 30) = 18.64
# 18.64 es el monto que debe pagar una persona de 52 años para recibir 30 pagos temporales anticipados 
# de un sol a partir de los 82 años.







# ANUALIDADES_TEMPORALES_VENCIDAS -----------------------------------------
#    ax:n| (con una rayita encima del n)

#¿Cuál es la relación entre ¨ax:n|y ax:n|? (con rayita encima de las n)

# ax:n| = ¨ax:n| + 1/(1 + i)^n × npx − 1

#Ejercicio:
#Una mujer sana de 64 años desea contratar un seguro de vida temporal de 27
#años a fin de recibir anualidades de 12 mil soles al final de cada año, considerando una
#TEA de 4%, calcule el valor de la prima bajo las siguientes condiciones

# Suponga que la mujer vivirá, con certeza, 91 años. ----------------------
i = 0.04
v = 1/(1+i)
d = i/(1+i)
n = 27
(a = (1 - v^n)/i)
a*12000

#  a27| = 16.33
# Suponga que la mujer está sujeta a morir según la tabla de mortalidad de la SBS --------

Tabla = probs2lifetable(probs = SPPS2017M*(1-AaxM)^4, radix = 10^6, type = "qx") 
(P = axn(Tabla, x=64, n=27, i=0.04, payment = "arrears"))
#[1] 14.18596

#otras formas: 

(P = axn(Tabla, x=64, n=27, i=0.04) + 1/1.04^27*pxt(Tabla,x=64,t=27) - 1)
# [1] 14.18596
(P = 1/0.04*(1-1.04*AExn(Tabla, x=64, n=28, i=0.04)))
# [1] 14.18596

# a64:27| = 14.19



















# ANUALIDADES_TEMPORALES_CONTINUAS ----------------------------------------
#   ax:n| (con rayita encima del a y n)

#Calcule el valor de la prima que debe pagar una persona de 43 años cuyo
#tiempo de vida puede ser explicado mediante la función de supervivencia:

#                S(x) = 1 − (x^2.1/17209.57)

# para 0 <= x <= 104, a fin de recibir 20 mil soles anuales durante 20 años. Considere una TEA de 2%.

S = function(x){1-(x/104)^2.1}
d = log(1+0.02)
f = function(t){exp(-d*t)*S(43+t)/S(43)}
(P = integrate(f,0,20)$value) # aquí 104-43 es reemplazado por 20
## [1] 14.88378
P*20000
## [1] 297675.7




# ANUALIDADES_DIFERIDAS ---------------------------------------------------

# m|¨ax = ANUALIDAD VITALICIA ANTICIPADA, diferida "m" años ------------

#                           m|¨ax = mEx * ¨ax+m
#                           m|¨ax = ¨ax − ¨ax:m|
#                     Diferi.Anti = V.anticipados - Temporales

#                          ¨ax = ¨ax:m| + m|¨ax

# Ejemplo:
# Un hombre inválido de 35 años contrata un seguro, mediante el cual recibirá
# beneficios de 15 mil soles al inicio de cada año, iniciando de aquí a 30 años y hasta que
# fallezca. Considerando una tasa anual de 1%, calcule el valor de la prima bajo las
# siguientes condiciones:

# Suponga que vivirá, con certeza, 110 años (serían 45 pagos)

i = 0.01
v = 1/(1+i)
d = i/(1+i)
n = 45
(a = (1 - v^n)/i)
## [1] 36.09451
a*15000

# Interpretacion: ¨a45|(rayita encima del 45) = 36.1
# 36.1 es el monto que hoy hay que pagar para recibir 45 pagos diferidos de 1 sol


# Suponga que está sujeto a morir según la tabla de mortalidad de la SBS.

Tabla = probs2lifetable(probs = SPPI2017H*(1-AaxH)^4, radix = 10^6, type = "qx")

#primera forma: m|¨ax = mEx * ¨ax+m

x = 65
i = 0.01

(a65 = axn(Tabla, x , i = i))

x = 35
n = 30
i = 0.01
(E35 = Exn(Tabla, x = x, n = n, i = i))

a65*E35


#segunda forma:  m|¨ax = ¨ax − ¨ax:m|

x = 35
i = 0.01
n = 30
(P = axn(Tabla, x=x, i=i) - axn(Tabla, x=x, n=n, i=0.01))

#tercera forma:
x=35
m=30
i=0.01

(P = axn(Tabla, x=x, m=m, i=i))

P*15000

#interpretación: 30|¨a35 = 6.72
# 6.72 es el monto que hoy hay que pagar para recibir 1 sol de pagos vitalicios diferidos
#anticipados por 30 años a partir de los 65 años.


# m|ax  = ANUALIDAD VITALICIA VENCIDA ---------------------------------

#Ejemplo: Una mujer sana de 39 años contrata hoy un seguro vitalicio mediante el cual
#recibirá beneficios de 16 mil soles anuales a partir del 2038 (dentro de 17 años).
#Considerando una tasa anual de 1.3%, calcule el valor de la prima bajo suponiendo que
#está sujeta a morir según la tabla de mortalidad de la SBS.


#16|a39 = 17|¨a39

#                               m|ax = mEx * ax+m

x = 39 #edad
n1 = 16
n2 = 17
i = 0.013
m = 17 #años diferidos 

Tabla = probs2lifetable(probs = SPPS2017M*(1-AaxM)^4, radix = 10^6, type = "qx")

(P = Exn(Tabla,x=x,n=n1,i=i) * axn(Tabla,x=(x + n1), i=i, payment="arrears"))
# 21.14091

(P = Exn(Tabla, x=x, n=n2, i=i) * axn(Tabla, x= (x + n2) ,i=i))
# 21.14091

(P = axn(Tabla, x = x, i = i) - axn(Tabla,x = x, n = n2, i = i))
# 21.14091

(P = axn(Tabla, x=x, m=m, i=i))
# 21.14091

# m|¨ax:n| (con rayita encima de la n) = ANUALIDAD TEMPORAL ANTICIPADA --------

#                              m|¨ax:n| = mEx * ¨ax+m:n|

#Ejercicio: 
#Un hombre sano de 34 años contrata un seguro temporal de 21 años, diferido
#en 23 años, mediante el cual recibirá beneficios de 12500 soles al inicio de cada
#si hoy 18 de enero del 2021 firma la póliza, el primer pago sería el 18 de enero del 2044 y
#el último el 18 de enero del 2064). Considerando una tasa anual de 2.6%, calcule el
#valor de la prima suponiendo que está sujeto a morir según las tablas de la SBS.

x = 34 #edad
n = 23 #número de anualidades
i = 0.026
m = 23 #años diferidos 

Tabla = probs2lifetable(probs = SPPS2017H*(1-AaxH)^4, radix = 10^6, type = "qx")

(P = axn(Tabla,x = x, n=n ,m=m, i=i, payment="arrears"))

P*12500

#  m|ax:n| (con rayita encima de la n)  = ANUALIDAD TEMPORAL VENCIDA --------

#                              m|ax:n| = mEx × ax+m:n|

# usar la función axn con los argumentos m, n y arrears












#                              m|ax:n| = mEx * ax+m:n|





# ANUALIDADES DIFERIDAS CONTINUAS -----------------------------------------

# 15|a38 (rayita encima de la a)

#Ejemplo:

#Calcule el valor de la prima que debe pagar una persona de 38 años cuyo
#tiempo de vida puede ser explicado mediante la función de supervivencia:

#                     S(x) = 1 − (x^2.1/17209.57)

#para 0 <= x <= 104, a fin de recibir 20 mil soles anuales, siendo el primer pago dentro de
#15 años. Considere una TEA de 3.9%.

S = function(x){1-(x/104)^2.1}
d = log(1+0.02)
f = function(t){exp(-d*t)*S(38+t)/S(38)}
(P = integrate(f,15,104-38)$value) # aquí 0 es reemplazado por 15
# 13.16993
P*20000
# 263398.6





# ANÁLISIS DE MULTIPLES TIEMPOS DE VIDA -----------------------------------

# PROBABILIDADES INMEDIATAS -----------------------------------------------

# SUPERVIVENCIA CONJUNTA --------------------------------------------------

#                        npxy = npx×npy

#Calcule la probabilidad de que en el 2021, una mujer sana de 40 años y un hombre
#inválido de 43 años vivan al menos 30 años más.

Tabla1 = probs2lifetable(probs = SPPS2017M*(1-AaxM)^4, 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 1")
Tabla2 = probs2lifetable(probs = SPPI2017H*(1-AaxH)^4, 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 2")
Tablas = list(T1 = Tabla1, T2 = Tabla2)
pxyzt(Tablas, x = c(40,43), t = 30, status = "joint")

#otra forma

pxt(Tabla1, x=40, t=30)*pxt(Tabla2, x=43, t=30)


# DISOLUCIÓN --------------------------------------------------------------

#                          nqxy = 1−npxy

#Una mujer sana tiene 40 años y un hombre inválido, 43 años. Calcule la probabilidad
#de que al menos uno de los asegurados fallezca antes de los próximos 30 años (a partir del 2021)

1-pxyzt(Tablas, x = c(40,43), t = 30, status = "joint")

#otra forma
qxyzt(Tablas, x = c(40,43), t = 30, status = "joint")


# NO EXTINCIÓN ------------------------------------------------------------

#                         npxy = npx+npy−npxy

#Una mujer sana tiene 40 años y un hombre inválido, 43 años. Calcule la probabilidad de que al menos
#uno de los asegurados viva al menos 30 años más.

npx  = pxt(Tabla1, x=40, t=30)
npy  = pxt(Tabla2, x=43, t=30)
npxy = pxyzt(Tablas, x=c(40,43), t=30, status = "joint")
npx + npy - npxy

#otra forma

npx + npy - npx*npy

#otra forma

pxyzt(Tablas, x=c(40,43), t=30, status = "last")


# DISOLUCIÓN Y NO EXTINCIÓN -----------------------------------------------

#                       nPxy(1)= nPx+nPy−2nPxy = nPxy − nPxy

#Una pareja está conformada por una mujer sana de 40 años y un hombre inválido de 43 años. 
#Calcule la probabilidad de exactamente uno de ellos sobreviva al menos 30 años más.

npx  = pxt(Tabla1, x = 40, t = 30)
npy  = pxt(Tabla2, x = 43, t = 30)
npxy = pxyzt(Tablas, x=c(40,43), t = 30, status = "joint")
npx + npy - 2*npxy

#otra forma

npx + npy - 2*npx*npy


# EXTINCIÓN ---------------------------------------------------------------

#                         nQxy¯¯¯¯¯¯= 1−nPxy¯¯¯¯¯¯

#Una pareja está conformada por una mujer sana de 40 años y un hombre inválido de 43 años.
#Calcule la probabilidad de que el último de ellos fallezca antes de los próximos 30 años.

npx  = pxt(Tabla1, x=40, t=30)
npy  = pxt(Tabla2, x=43, t=30)
npxy = pxyzt(Tablas, x=c(40,43), t=30, status = "joint")
1-(npx + npy - npxy)

#otra forma
1-(npx + npy - npx*npy)

#otra forma
qxyzt(Tablas, x=c(40,43), t=30, status = "last")












# PROBABILIDADES DIFERIDAS ------------------------------------------------

#Se asume que transurrirán "m" años, luego de eso puede pasar dos cosas
# 1) Ambas personas siguen vivas, con probabilidad mPx,y
# 2) Al menos unas de las personas falleció (solo vive uno ninguno) con proba = mQx,y


# DISOLUCIÓN --------------------------------------------------------------

#                m|nQxy = mPxy * nQx+m;y+m

#Una pareja de esposos tiene 40 (esposa inválida) y 47 (esposo sano) años respectivamente. 
#Calcule la probabilidad de que ambos asegurados vivan al menos 28 años a partir del 2021, 
#pero luego al menos uno de ellos fallezca antes de los 20 siguientes:

Tabla3 = probs2lifetable(probs = SPPI2017M*(1-AaxM)^4, 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 3")
Tabla4 = probs2lifetable(probs = SPPS2017H*(1-AaxH)^4, 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 4")
Tablas2 = list(T3 = Tabla3, T4 = Tabla4)

(p1  = pxyzt(Tablas2, x=c(40,47), t=28, status = "joint"))

(q2  = qxyzt(Tablas2, x=c(40+28,47+28), t=20, status = "joint"))

p1*q2

#es la probabilidad de que ambos vivan al menos 28 años más pero luego de dicho periodo al menos
#uno de ellos fallezca antes de los próximos 20 años.


# EXTINCIÓN ---------------------------------------------------------------

#                      m|nQxy¯¯¯¯¯¯ = mPxy * nQx+m;y+m¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯

#Una pareja de esposos tiene 40 (esposa inválida) y 47 (esposo sano) años respectivamente. 
#Calcule la probabilidad de que ambos asegurados vivan al menos 28 años a partir del 2021, 
#pero luego el último de ellos fallezca antes de los 20 siguientes:

(p1  = pxyzt(Tablas2, x=c(40,47), t=28, status = "joint"))
(q2  = qxyzt(Tablas2, x=c(40+28,47+28), t=20, status = "last"))
p1*q2

#es la probabilidad de que ambos vivan al menos 28 años más pero luego de dicho periodo ambos 
#fallezcan (o el último sobreviviente fallezca) antes de los próximos 20 años.











# SEGUROS -----------------------------------------------------------------

# SEGURO VITALICIO CONJUNTO -----------------------------------------------

#Cuando fallezca uno, deja de ser conjunta.
#la indemnización se paga cuando fallece el primer asegurado 
#el valor esperado actuarial del seguro de vida completa es:

#                   Axy = E(Zxy) = v^k * k|Qxy

#Una pareja de esposos sanos de 50 (esposo) y 53 años (esposa) contrata un seguro de vida conjunta completa. 
#Determine la prima a ser pagada considerando una tasa de interés de 2% anual, para una cobertura de 200 mil soles.

Tabla5 = probs2lifetable(probs = SPPS2017M*(1-AaxM)^4, 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 5")
Tabla6 = probs2lifetable(probs = SPPS2017H*(1-AaxH)^4, 
                         radix = 10^6,   
                         type = "qx",  
                         name = "Tabla 6")
Tablas3 = list(T5 = Tabla5, T6 = Tabla6)
A = Axyzn(Tablas3, x = c(53,50), i = 0.02, status = "joint")
A  #esta es la prima unitaria 

S = 200000
P = S*A
P #monto de la prima única


# SEGURO VITALICIO DEL ULTIMO SOBREVIVIENTE -------------------------------

#la indemnizacion se paga cuando han fallecido ambos asegurados
#el valor esperado acturial del segundo de vida completa es:

#        Axy¯¯¯¯¯¯ = E(Zxy¯¯¯¯¯¯) = ∑ ν^k * k|Qxy¯¯¯¯¯¯


#Una pareja de esposos sanos de 50 (esposo) y 53 años (esposa) contrata un seguro de vida completa para el último sobreviviente.
#Determine la prima a ser pagada considerando una tasa de interés de 2% anual, para una cobertura de 200 mil soles.

A = Axyzn(Tablas3, x = c(53,50), i = 0.02, status = "last")
A #prima unitaria 
S = 200000
P = S*A
P # prima unica a pagar 

































