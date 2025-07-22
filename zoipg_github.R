########### R code for data analysis##########
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



###### R code for zoipg simulation #############
library(rootSolve)
sim=1000
n=100
p0=0.3
p1=0.9
theta=0.5
p=(theta+1)/(theta+2)
MM1=MM2=MM3=c();
ML1=ML2=ML3=c();
EM1=EM2=EM3=c();
for(i in 1:sim)
 {
x1=rexp(n,theta)
x2=rgamma(n,rate=theta,shape=2)
mix=runif(n)
data=ifelse(mix<p,x1,x2)
V=rpois(n,data)
  B1=rbinom(n,1,p0)
  B2=rbinom(n,1,p1)
  y=B1*(1-B2)+(1-B1)*V
mean(y)
myu=p0*(1-p1)+(1-p0)*(theta+3)/(theta*(theta+2))

# MOMENT ESTIMATION
fm1=mean(y)
  fm2=mean(y*(y-1))
  fm3=mean(y*(y-1)*(y-2))
s1=fm3/fm2
MM3[i]=((3-4*s1)+sqrt((4*s1+3)^2+12*s1))/(2*s1)
  MM1[i]=1-(fm2*(MM3[i]+2)*MM3[i]^2)/(2*(MM3[i]+4))
 muhat=(MM3[i]+3)/(MM3[i]*(MM3[i]+2))
p21=MM1[i]-fm1+(1-MM1[i])*muhat
  p22=p21/MM1[i]
  MM2[i]=p22;
if ( MM1[i]< 0 || MM1[i]> 1  )
{MM1[i]=0.5
}
if( MM2[i]<0 || MM2[i]>1)
{MM2[i]=0.1}

### Likelihood estimation ##########
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
 ML3[i]=uniroot(loglt,c(0,100))$root
loglp=function(x)
     {
      p00=x[1]
      p10=x[2]
      c1=(ML3[i]*(ML3[i]^2+3*ML3[i]+1))/((ML3[i]+2)*(ML3[i]+1)^2)
      c2=(ML3[i]*(ML3[i]^2+4*ML3[i]+1))/((ML3[i]+2)*(ML3[i]+1)^3)

      F1=c1*(1-p00)+p00*p10-s0/n
      F2=c2*(1-p00)+p00*(1-p10)-s1/n
      return(c(F1=F1,F2=F2))
     }
  ML1[i]=multiroot(loglp,c(0.5,0.5))$root[1]
  ML2[i]=multiroot(loglp,c(0.5,0.5))$root[2]


p01=(s0+s1)/n
p02num= ML3[i]*( ML3[i]^3+5* ML3[i]^2+8* ML3[i]+2)
p02denom=( ML3[i]+2)*( ML3[i]+1)^3
ML1[i]=p01-(p02num/p02denom)
p11=s0/n
p12=(1- ML1[i])*( ML3[i]/( ML3[i]+2))*(( ML3[i]^2+3* ML3[i]+1)/(( ML3[i]+1)^2))
ML2[i]=(p11-p12)/ ML1[i]
if ( ML1[i]< 0 || ML1[i]> 1  )
{ML1[i]=0.5
}
if( ML2[i]<0 || ML2[i]>1)
{ML2[i]=0.1}

#####EM algorithm ####
 EMfit=function(par)
  {
   EM=function(para)
      {
	 EMp0=para[1];
	 EMp1=para[2];
	 EMtheta=para[3];
	 a1=(EMp0*EMp1)/(EMp0*EMp1+(1-EMp0)*(EMtheta*(EMtheta^2+3*EMtheta+1))/((EMtheta+2)*(EMtheta+1)^2))
	 a2=(EMp0*(1-EMp1))/(EMp0*(1-EMp1)+(1-EMp0)*(EMtheta*(EMtheta^2+4*EMtheta+1))/((EMtheta+2)*(EMtheta+1)^3))
	 B1em=a1*(y==0)+a2*(y==1)+0*(y>=2)
	 b11=EMp0*EMp1+((1-EMp0)*EMp1*EMtheta*(EMtheta^2+3*EMtheta+1))/((EMtheta+2)*(EMtheta+1)^2)
	 b12=EMp0*EMp1+((1-EMp0)*(EMtheta*(EMtheta^2+3*EMtheta+1)))/((EMtheta+2)*(EMtheta+1)^2)
	 b1=b11/b12
	 b21=((1-EMp0)*EMp1*EMtheta*(EMtheta^2+4*EMtheta+1))/((EMtheta+2)*(EMtheta+1)^3)
	 b22=EMp0*(1-EMp1)+((1-EMp0)*EMtheta*(EMtheta^2+4*EMtheta+1))/((EMtheta+2)*(EMtheta+1)^3)
	 b2=b21/b22
	 b3=EMp1
	 B2em=b1*(y==0)+b2*(y==1)+b3*(y>=2)
	  logtheta=function(th)
	   {
	    -sum((1-B1em)*log(th*(th*y+th^2+3*th+1)/((th+2)*(th+1)^(y+2))))
	    }
	 thetaEM=optimize(logtheta,c(0,100))$minimum
	 para[1]=mean(B1em)
	 para[2]=mean(B2em)
	 para[3]=thetaEM
       
       loglik=function(theta)
	  {
         p01=theta[1]
         p11=theta[2]
         th=theta[3]
if(p01<0 || p01>1)
{p01=0.5}
if(p11<0 || p11>1)
{p11=0.5}
         logL= sum(B1em*log(p01)+(1-B1em)*(log(1-p01)+2*log(th)+
               log(y+th+2)-(y+3)*log(th+1))+B2em*log(p11)+
               (1-B2em)*log(1-p11))
         return(logL)
         }
   list(pa=para,ll=loglik(para))
 }
iter=0;
para.old=par;
like.old=EM(par)$ll;
repeat
 {
  iter=iter+1
  para.new <- EM(para.old)$pa
  loglik.new <- EM(para.new)$ll
  if(abs(loglik.new-like.old)<0.001)
  break
  else{para.old=para.new
       like.old=loglik.new
       }
    }
   list(para=para.new,loglik=loglik.new)
  }
 mm=c(MM1[i],MM2[i],MM3[i])
EM1[i]=EMfit(mm)$para[1]
EM2[i]=EMfit(mm)$para[2]
EM3[i]=EMfit(mm)$para[3]
}

sum(EM2>0)
True=c(theta,p0,p1)
MM=c(mean(MM3),mean(MM1),mean(MM2))
MLE=c(mean(ML3),mean(ML1),mean(ML2))
EM=c(mean(EM3),mean(EM1),mean(EM2))
MSEmm=c(mean((MM3-theta)^2),mean((MM1-p0)^2),mean((MM2-p1)^2))
MSEml=c(mean((ML3-theta)^2),mean((ML1-p0)^2),mean((ML2-p1)^2))
MSEem=c(mean((EM3-theta)^2),mean((EM1-p0)^2),mean((EM2-p1)^2))
cbind(True,MM,MSEmm,MLE,MSEml,EM,MSEem)
n
p0
p1
theta





