# Temas Selectos de Probabilidad y Estadística

# Segundo Examen Parcial

### 16 de marzo
Analizar la deducción de la cópula subyacente en el Ejemplo 3.3 del libro de Nelsen y calcular en función del parámetro de dicha cópula las medidas de dependencia de Schweizer-Wolff, Hoeffding y distancia supremo, así como las medidas de concordancia de Kendall, Spearman y Erdely. Luego graficar juntas todas las medidas de dependencia y concordancia como funciones del parámetro.

```R
#Cargar librería "cubature" para integrar
install.packages("cubature")
library(cubature)
```

Medida de dependencia Schweizer-Wolff:
```R
Medida1<-function(a)
{
  SW<-function(x) 
    {
      if(0<=x[1] && x<=a*x[2] && a*x[2]<=a) C<-x[1]
      if(0<=a*x[2] && a*x[2]<x[1] && x[1]<1-(1-a)*x[2]) C<-a*x[2]
      if(a<=1-(1-a)*x[2] && 1-(1-a)*x[2]<=x[1] && x[1]<=1) C<-x[1]+x[2]-1
      abs(C-x[1]*x[2])
    }

  12*adaptIntegrate(SW, lowerLimit = c(0, 0), upperLimit = c(1, 1),maxEval=500)$integral
}

Medida1a<-Vectorize(Medida1)
curve(Medida1a,0,1,xname="SW",xlab="Medida de dependencia Schweizer-Wolff")
```
![SW](images/Medida1.png)

Medida de dependencia Hoeffding:
```R

Medida2<-function(a)
{ 
  Hoeffding<-function(x) 
    {
      if(0<=x[1] && x[1]<=a*x[2] && a*x[2]<=a) C<-x[1]
      if(0<=a*x[2] && a*x[2]<x[1] && x[1]<1-(1-a)*x[2]) C<-a*x[2]
      if(a<=1-(1-a)*x[2] && 1-(1-a)*x[2]<=x[1] && x[1]<=1) C<-x[1]+x[2]-1
      (C-x[1]*x[2])^2
    }
  
  sqrt(90*adaptIntegrate(Hoeffding, lowerLimit = c(0, 0), upperLimit = c(1, 1),maxEval=500)$integral)
}

Medida2a<-Vectorize(Medida2)
curve(Medida2a,0,1,xname="H",xlab="Medida de dependencia Hoeffding")
```
![H](images/Medida2.png)

Medida de dependencia distancia supremo:
```R
Medida3<-function(a)
{ 
  Supremo<-function(x) 
    {
      if(0<=x[1] && x[1]<=a*x[2] && a*x[2]<=a) C<-x[1]
      if(0<=a*x[2] && a*x[2]<x[1] && x[1]<1-(1-a)*x[2]) C<-a*x[2]
      if(a<=1-(1-a)*x[2] && 1-(1-a)*x[2]<=x[1] && x[1]<=1) C<-x[1]+x[2]-1
      abs(C-x[1]*x[2])
    }
  
  4*optim(c(.5,.5),Supremo,control=list(fnscale=-1))$value
  }

Medida3a<-Vectorize(Medida3)
curve(Medida3a,0,1,n=300,xname="Sup",xlab="Medida de dependencia Supremo")
```
![Sup](images/Medida3.png)

Medida de concordancia de Kendall:
```R
#Concordancia1
Concordancia1<-function(a)
{
  Kendall<-function(x) 
  {
    1*(a<=1-(1-a)*x[2])*(1-(1-a)*x[2]<=x[1])*(x[1]<=1)
  }
  1-4*adaptIntegrate(Kendall, lowerLimit = c(0, 0), upperLimit = c(1, 1),maxEval=1000)$integral
}

Concordancia1a<-Vectorize(Concordancia1)
curve(Concordancia1a,0,1,xname="K",xlab="Medida de concordancia Kendall")
```
![Ke](images/Concordancia1.png)

Medida de concordancia de Spearman:
```R
Concordancia2<-function(a)
{
  Sp<-function(x) 
  {
    if(0<=x[1] && x<=a*x[2] && a*x[2]<=a) C<-x[1]
    if(0<=a*x[2] && a*x[2]<x[1] && x[1]<1-(1-a)*x[2]) C<-a*x[2]
    if(a<=1-(1-a)*x[2] && 1-(1-a)*x[2]<=x[1] && x[1]<=1) C<-x[1]+x[2]-1
    C-x[1]*x[2]
  }

  12*adaptIntegrate(Sp, lowerLimit = c(0, 0), upperLimit = c(1, 1),maxEval=500)$integral
}

Concordancia2a<-Vectorize(Concordancia2)
curve(Concordancia2a,0,1,xname="Sp",xlab="Medida de concordancia Spearman")
```
![Sp](images/Concordancia2.png)

Medida de concordancia de Erdely:
```R
Concordancia3<-function(a)
{ 
  Erdely1<-function(x) 
    {
      1*(a<=1-(1-a)*x[2])*(1-(1-a)*x[2]<=x[1])*(x[1]<=1)-x[1]*x[2]
    }
  
  Erdely2<-function(x) 
    {
      x[1]*x[2]-1*(a<=1-(1-a)*x[2])*(1-(1-a)*x[2]<=x[1])*(x[1]<=1)
    }
  
  4*(optim(c(a,a),Erdely1,control=list(fnscale=-1))$value
     - optim(c(1-a,1-a),Erdely2,control=list(fnscale=-1))$value)
}

Concordancia3a<-Vectorize(Concordancia3)
curve(Concordancia3a,0,1,xname="Er",xlab="Medida de concordancia Erdely")
```
![Er](images/Concordancia3.png)


### 17 de marzo
Considere un vector aleatorio (X,Y) con función de densidad conjunta de probabilidades del Ejemplo 1.7 de las notas sobre vectores aleatorios. Programando en R:
a) Simule una muestra aleatoria de (X,Y) de tamaño n = 3000 y realice un gráfico de dispersión.
```R
ui<-runif(3000,0,1)
xi<--log(1-ui)
vi<-runif(3000,0,1)
yi<-xi-log(1-vi)
plot(xi,yi)
```
b) Con los valores simulados de X obtenga un histograma en la escala adecuada para que encima grafique la densidad teórica marginal de X. Lo mismo para Y.
c) Obtenga gráficas de los conjuntos de nivel de la cópula subyacente C mediante las funciones contour e image.
d) Igual que en el inciso anterior pero de C(u,v) - uv.
e) Calcule las medidas de dependencia de Schweizer-Wolff, Hoeffding y distancia supremo, así como las medidas de concordancia de Kendall, Spearman y Erdely.
```R
#b) Con los valores simulados de X obtenga un histograma en la 
#escala adecuada para que encima grafique la densidad teórica 
#marginal de X. Lo mismo para Y.
hist(xi,prob=T)
curve(dexp(x),from=min(xi),to=max(xi),add=T)

?uniroot()

hist(yi,prob=T)
curve(dgamma(x,2,1),from=min(yi),to=max(yi),add=T)

```
![Sexo](images/sexo.png)

### 25 de marzo
Simule una muestra aleatoria de tamaño n = 3000 a partir de un vector aleatorio (U,V) con marginales Uniformes(0,1) y cópula Clayton con parámetro 2. Realice gráficas de histogramas marginales y uno de dispersión. Calcule las medidas de dependencia de Schweizer-Wolff, Hoeffding y distancia supremo, así como las medidas de concordancia de Kendall, Spearman y Erdely.
```R
ui<-runif(10)->xi
vi<-runif(10)
#vi<-
  
Clayton<-function(u,v)
{
  (max(u^(-2)+v^(-2)-1,0))^(-1/2)
}

FVlU<-function(v,u)
{
  ((max(u^(-2)+v^(-2)-1,0))^(-1/2))/u
}

FV_inv <- function(v,u){
  
  uniroot(function(x) FVlU(x,u)-v, interval = c(0,1))$root
  
}

yi<-vector(mode='numeric',length = length(ui))

for(i in 1:length(ui)) yi[i]<-FV_inv(vi[i],xi[i])

plot(ui,yi)

concor<-sum(outer(yi,yi,function(x,y) x-y)*outer(ui,ui,function(x,y) x-y)>0)
discor<-sum(outer(yi,yi,function(x,y) x-y)*outer(ui,ui,function(x,y) x-y)<0)

(concor-discor)/(concor+discor)

```
![tabla](images/edad.png)

### 25 de marzo
Lo mismo del ejercicio anterior pero para (X,Y) pero con distribuciones marginales X ~ Pareto(1,1), Y ~ Cauchy(0,1). Calcule las medidas de dependencia de Schweizer-Wolff, Hoeffding y distancia supremo, así como las medidas de concordancia de Kendall, Spearman y Erdely.
```R
#Código
```
![tabla](images/edad.png)
       
