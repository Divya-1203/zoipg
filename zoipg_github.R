library(rootSolve)
x=c(0,1,2,3,4,5,6,7,8,9,14)
n1=length(x)
o=c(64,55,22,12,6,3,2,1,1,1,1)
y=c(rep(x,o))
n=length(y)
z=table(y)
z
o<- as.numeric(z)
x <- as.numeric(names(z))
n1=length(x)
zero_proportion=sum(y==0)/n
k1=k2=k3=E=e=p=c()
mean(y)
var(y)
dis=var(y)/mean(y)
dis
ep=o/n
ep
########ZERO-ONE INFLATED POISSON####################################################

###########maximum likelihood estimation######

s0=sum(y==0)
s1=sum(y==1)
s2=n-s0-s1
y1=y[y>=2]
su= sum(y1)

loglt=function(th)
{
F=su*(exp(th)-th-1)-s2*th*(exp(th)-1)
return(F=F) 
}
ML3=uniroot(loglt,lower = 0.0005,upper= 100)$root
q0=s0/n
q1=s1/n
d1=1-exp(-ML3)-ML3*exp(-ML3)
d1

ML1=(q0+q1-(exp(-ML3)*(1+ML3)))/(d1)
ML2= (q0-(1-ML1)*exp(-ML3))/(ML1)
ML1
ML2
ML3

#######AIC#############
m=3
 ###########log-likelihood###############

u1=s0*log(ML1*ML2+(1-ML1)*exp(-ML3))
u2=s1*(log(ML1*(1-ML2)+(1-ML1)*ML3*exp(-ML3)))
u3=sum(log((1-ML1)*exp(-ML3)*(ML3)^y1/factorial(y1)))
LL=u1+u2+u3
AICML=2*m-2*LL
AICML
 ##############BIC###############
BICML= m*log(n)-2*LL
BICML
mmean=ML1*(1-ML2)+(1-ML1)*ML3
mmean
mm2=ML1*(1-ML2)+(1-ML1)*(ML3+ML3^2)
vvar=mm2-mmean^2
vvar
disp=vvar/mmean
disp


####ZERO-ONE INFLATED POISSON GARIMA########################################

##########moment estimation########

fm1=mean(y)
  fm2=mean(y*(y-1))
  fm3=mean(y*(y-1)*(y-2))
s1=fm3/fm2
MM3=((3-4*s1)+sqrt((4*s1+3)^2+12*s1))/(2*s1)
  MM1=1-(fm2*(MM3+2)*MM3^2)/(2*(MM3+4))
 muhat=(MM3+3)/(MM3*(MM3+2))
p21=MM1-fm1+(1-MM1)*muhat
  p22=p21/MM1
  MM2=p22;


###########maximum likelihood estimation#####
 s0=sum(y==0)
  s1=sum(y==1)
  y1=y[y>=2]
  loglt=function(x)
      {
       a1=n-s0-s1
       a2=(4*x^2+4*x+2)/(x*(x+1)*(x^2+5*x+2))
       a3=sum((2*x+3+y1)/(x^2+3*x+1+x*y1))
       a4=sum((y1)/(x+1))
       F=a1*a2+a3-a4;
       return(F=F)  
      }
 ML3=uniroot(loglt,c(0,100))$root
loglp=function(x)
     {
      p00=x[1]
      p10=x[2]
      c1=(ML3*(ML3^2+3*ML3+1))/((ML3+2)*(ML3+1)^2)
      c2=(ML3*(ML3^2+4*ML3+1))/((ML3+2)*(ML3+1)^3)

      F1=c1*(1-p00)+p00*p10-s0/n
      F2=c2*(1-p00)+p00*(1-p10)-s1/n
      return(c(F1=F1,F2=F2))
     }
  ML1=multiroot(loglp,c(0.5,0.5))$root[1]
  ML2=multiroot(loglp,c(0.5,0.5))$root[2]
MM1
MM2
MM3
ML1
ML2
ML3

#######AIC#############
m=3
 ###########log-likelihood###############
u1=s0*log(ML1*ML2+(1-ML1)*(ML3*(ML3^2+3*ML3+1))/((ML3+2)*(ML3+1)^2))
u2=s1*log(ML1*(1-ML2)+(1-ML1)*(ML3*(ML3^2+4*ML3+1))/((ML3+2)*(ML3+1)^3))
u3=(n-s0-s1)*log((ML3*(1-ML1))/((ML3+2)*(ML3+1)^2))
u4=sum(log(ML3^2+3*ML3+1+ML3*y1)-y1*log(ML3+1))
LL=u1+u2+u3+u4
AICML=2*m-2*LL
AICML

############BIC#####################
BICML= m*log(n)-2*LL
BICML
mmean=ML1*(1-ML2)+(1-ML1)*((ML3+3)/(ML3*(ML3+2)))
mmean
mm2=ML1*(1-ML2)+(1-ML1)*((ML3^2+5*ML3+8)/(ML3^2*(ML3+2)))
vvar=mm2-mmean^2
vvar
disp=vvar/mmean
disp

###########ZERO ONE INFLATED POISSON LINDLEY#######################################

### Likelihood estimation ##########
s0=sum(y==0)
s1=sum(y==1)
y1=y[y>=2]
loglt=function(x){
       a1=n-s0-s1
       a2=(x^3+8*x^2+7*x+2)/(x*(x+1)*(x^2+4*x+1))
       a3=sum(1/(x+2+y1))
       a4=1/(x+1)*sum(y1)
       F=a1*a2+a3-a4;
       return(F=F)  
}
ML3=uniroot(loglt,c(0,100))$root
loglp=function(x){
p0=x[1]
p1=x[2]
c1=ML3^2*(ML3+2)/(ML3+1)^3
c2=ML3^2*(ML3+3)/(ML3+1)^4
F1=c1*(1-p0)+p0*p1-s0/n
F2=c2*(1-p0)+p0*(1-p1)-s1/n
return(c(F1=F1,F2=F2))
}
library(rootSolve)
ML1=multiroot(loglp,c(0.5,0.5))$root[1]
ML2=multiroot(loglp,c(0.5,0.5))$root[2]
ML1
ML2
ML3

#######AIC#############
m=3
u1=s0*log(ML1*ML2+(1-ML1)*(ML3^2*(ML3+2))/((ML3+1)^3))
u2=s1*log(ML1*(1-ML2)+(1-ML1)*(ML3^2*(ML3+3))/((ML3+1)^4))
u3=(n-s0-s1)*log((1-ML1)*(ML3^2)/(ML3+1)^3)
u4=sum(log(ML3+2+y1)-y1*log(ML3+1))
LL=u1+u2+u3+u4
AICML=2*m-2*LL
AICML
###########BIC##############
BICML= m*log(n)-2*LL
BICML
mmean=ML1*(1-ML2)+(1-ML1)*((ML3+2)/(ML3*(ML3+1)))
mmean
mm2=ML1*(1-ML2)+(1-ML1)*((ML3^2+4*ML3+6)/(ML3^2*(ML3+1)))
vvar=mm2-mmean^2
vvar
disp=vvar/mmean
disp



#

