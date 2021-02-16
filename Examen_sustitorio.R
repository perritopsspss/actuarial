
# EXAMEN SUSTITUTOTIO ACTUARIAL 


library(lifecontingencies)
TablaSBS = read.table("TMSPP2017.txt",T)

# Pregunta 1 --------------------------------------------------------------

#mujer invalida
probas = as.vector(TablaSBS$SPPI2017M*(1-TablaSBS$AaxM)^4)
#px= asociado a probabilidades de supervivencia, qx: asociado a probabildiades de fallecimiento
Tabla4 = probs2lifetable(probs = probas,
                         radix = 10^6,
                         type = "qx",
                         name = "Tabla4")
#17q31 = probabilidad de que una mujer invalida de 31 años viva por lo menos 17 años mas
qxt(Tabla4, 31, 17)

#m|nQx = mPx * nQx+m (probabilidad de que una mujer invalida de 31 años, alcance la edad 62 pero fallezca antes del siguiente)

#32|1q31 
x= 31
m=32
n=1

pxt(Tabla4, x = x, t = m) * qxt(Tabla4, x = x+m, t = n)
# Pregunta 2 --------------------------------------------------------------


# Pregunta 3 --------------------------------------------------------------


# Pregunta 4  -------------------------------------------------------------
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

# Pregunta 5  -------------------------------------------------------------

mu = function(s) return(0.3+0*s)
(I = integrate(mu,0,5)$value)
exp(-I)
# Pregunta 6 --------------------------------------------------------------


# Pregunta 7 --------------------------------------------------------------
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


# Pregunta 8  -------------------------------------------------------------

edad = 10.12
Sx = function(x) (1 - (x^0.25)/(105^0.25) )^0.1
f = function(m) return(Sx(edad+m)/Sx(edad)-0.5)
mediana = uniroot(f, lower = 0, upper = 105-edad)$root


num = integrate(Sx, lower = edad, upper = 105)$value
den = Sx(edad)
(e30 = num/den)

# Pregunta 9  -------------------------------------------------------------

#lx = l0(S(x))
100000*0.85


# Pregunta 10 -------------------------------------------------------------

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

q17-q18 #sale positivo.
