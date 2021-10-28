##Control parameters
#N=5: exponent on correlation is 1.2
#     multiple for y is 10
#N=10: exponent on correlation is 1.3
#     multiple for y is 10
#N=15: exponent on correlation is 1.3
#     multiple for y is 5
#N=20: exponent on correlation is 1.3
#     multiple for y is 5





start_time=Sys.time()
L=100;
par(mfrow=c(1,2))
M=20000;
burn=10000;
rho1=matrix(0,nrow=(M-burn),ncol=L); rho2=rho1;
G1=rho1; G2=rho1;
P1=rep(0,L); P2=P1;
max.d1=rep(0,L); max.d2=max.d1;
N=25;
x=seq(0,1,length.out=N);
q=function(x) {
  if(x==1) {
    return(.6);
  } else {
    return(.4);
  }
}
th1=function(x) {return(.25)}
th2=function(x) {exp(-(4*x-1)/4)}
eta=function(x,y1,y2) {
  result=y1*sin(2*pi*x)+y2*cos(2*pi*x);
  #result=y1*y2;
  return(result);
}


for (l in 1:L) {
  
  
  
  y=rnorm(N,mean=eta(x,th1(x),th2(x)),sd=0);
  
  R=function(p,x,x.p) {
    n=length(x);
    r=matrix(0,ncol=n,nrow=n);
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        r[i,j]=p^(abs(x[i]-x[j])^1);
      }
    }
    result=diag(n)+r+t(r);
    result=result
    return(result);
  }
  
  ####z[1:2N] theta, z[2N+1] rho1, z[2N+2] lambda1, z[2N+3] g1, z[2N+4] rho2, z[2N+5] lam2, z[2N+6] g2, z[2N+7] sig2
  
  e=.5; w=.1; alpha=2; beta=2; alpha1=2; beta1=.2;
  
  
  R.prop=R(.8,x,x); E.prop=eigen(R.prop); R.prop.5=E.prop$vectors%*%diag(sqrt(abs(E.prop$values)));
  th1.beg=as.vector(E.prop$vectors%*%diag(sqrt(abs(E.prop$values)))%*%as.matrix(rnorm(N,mean=0,sd=sqrt(c0[1]))))+.5;
  th2.beg=as.vector(E.prop$vectors%*%diag(sqrt(abs(E.prop$values)))%*%as.matrix(rnorm(N,mean=0,sd=sqrt(c0[3]))))+.5;
  while((max(th1.beg)>1) || (min(th1.beg)<0)) {
    th1.beg=as.vector(E.prop$vectors%*%diag(sqrt(abs(E.prop$values)))%*%as.matrix(rnorm(N,mean=0,sd=sqrt(c0[1]))))+.5;
  }
  while((max(th2.beg)>1) || (min(th2.beg)<0)) {
    th2.beg=as.vector(E.prop$vectors%*%diag(sqrt(abs(E.prop$values)))%*%as.matrix(rnorm(N,mean=0,sd=sqrt(c0[3]))))+.5;
  }
  c0=c(.01,.5,.01,.5,.1,.95,.9);
  
  
  acc=rep(0,7);
  z=c(th1.beg,th2.beg,.9,0,0,.9,0,0,0);
  MC.tot=matrix(0,ncol=2*N+7,nrow=M);
  pick.th=rbinom(1,1,.5);
  A1=R(z[2*N+1],x,x); E1.old=eigen(A1);
  A2=R(z[2*N+4],x,x); E2.old=eigen(A2);
  for ( i in 1:M) {
    zn=rnorm(2*N+3,mean=0,sd=sqrt(c(rep(c0[1],N),rep(c0[3],N),c0[2],c0[4],c0[5])));
    
    temp=z;
    if(((z[2*N+3]+z[2*N+6])==0) || ((z[2*N+3]+z[2*N+6])==2)) {
      asd=sample(2,1);
      if (asd==1) {
        temp[2*N+3]=1-z[2*N+3];
      } else {
        temp[2*N+6]=1-z[2*N+6];
      }
    } else {
      asd=sample(3,1,prob=c(.45,.45,.1));
      if (asd==1) {
        temp[2*N+3]=1-z[2*N+3];
      } else if (asd==2) {
        temp[2*N+6]=1-z[2*N+6];
      } else {
        temp[2*N+3]=1-z[2*N+3];
        temp[2*N+6]=1-z[2*N+6];
      }
    }
    
    
    ##choose which to update
    for (type in 1:2) {
      if (pick.th==1) { ###############################for theta 1
        
        temp[1:N]=as.vector(R.prop.5%*%as.matrix(zn[1:N]))+z[1:N];
        #while((max(temp[1:N])>1) || (min(temp[1:N])<0)) {
        #  temp[1:N]=as.vector(R.prop.5%*%as.matrix(rnorm(N,0,sqrt(c0[1]))))+z[1:N];
        #}
        
        u=log(runif(1,0,1));
        alph=logf1(4,temp,E1.old)-logf1(4,z,E1.old);
        #print(alph)
        if (alph>=u) {
          z[1:N]=temp[1:N];
          acc[1]=acc[1]+1;
        }
        ####add delete swap for rho and gamma
        a=c0[6];
        if( temp[2*N+3]==0) {
          temp[2*N+1]=runif(1,a,1);
        } else {
          temp[2*N+1]=runif(1,0,1);
        }
        A1=R(temp[2*N+1],x,x); E1.new=eigen(A1);
        u=log(runif(1,0,1));
        alph=logf1(1,temp,E1.new)-log(1/b*temp[2*N+3]+1/(1-a)*(1-temp[2*N+3]))-logf1(1,z,E1.old)+log(1/b*z[2*N+3]+1/(1-a)*(1-z[2*N+3]));
        if (alph>=u) {
          z[c(2*N+1,2*N+3)]=temp[c(2*N+1,2*N+3)];
          acc[6]=acc[6]+1;
          E1.old=E1.new;
        }
        
        ####M-H steps for variance parameters
        temp[2*N+2]=z[2*N+2]+zn[2*N+1];
        
        u=log(runif(1,0,1));
        alph=logf(3,temp,E1.old)-logf(3,z,E1.old)-temp[2*N+2]+z[2*N+2];
        if (alph>=u) {
          z[2*N+2]=temp[2*N+2];
          acc[2]=acc[2]+1;
        }
        
        pick.th=0; #z[2*N+2]=log(c);
        
        
      } else {  #############################for theta 2
        
        temp[(N+1):(2*N)]=as.vector(R.prop.5%*%as.matrix(zn[(N+1):(2*N)]))+z[(N+1):(2*N)];
        #while((max(temp[(N+1):(2*N)])>1) || (min(temp[(N+1):(2*N)])<0)) {
        #  temp[(N+1):(2*N)]=as.vector(R.prop.5%*%as.matrix(rnorm(N,0,sqrt(c0[3]))))+z[(N+1):(2*N)];
        #}
        u=log(runif(1,0,1));
        alph=logf1(7,temp,E2.old)-logf1(7,z,E2.old);
        #print(alph)
        if (alph>=u) {
          z[(N+1):(2*N)]=temp[(N+1):(2*N)];
          acc[3]=acc[3]+1;
        }
        ####add delete swap for rho and gamma
        a=c0[7];
        if( temp[2*N+6]==0) {
          temp[2*N+4]=runif(1,a,1);
        } else {
          temp[2*N+4]=runif(1,0,1);
        }
        A2=R(temp[2*N+4],x,x); E2.new=eigen(A2);
        u=log(runif(1,0,1));
        alph=logf1(5,temp,E2.new)-log(1/b*temp[2*N+6]+1/(1-a)*(1-temp[2*N+6]))-logf1(5,z,E2.old)+log(1/b*z[2*N+6]+1/(1-a)*(1-z[2*N+6]));
        if (alph>=u) {
          z[c(2*N+4,2*N+6)]=temp[c(2*N+4,2*N+6)];
          acc[7]=acc[7]+1;
          E2.old=E2.new;
        }
        
        ####M-H steps for variance parameters
        temp[2*N+5]=z[2*N+5]+zn[2*N+2];
        
        u=log(runif(1,0,1));
        alph=logf(6,temp,E2.old)-logf(6,z,E2.old)-temp[2*N+5]+z[2*N+5];
        if (alph>=u) {
          z[2*N+5]=temp[2*N+5];
          acc[4]=acc[4]+1;
        }
        
        pick.th=1; #z[2*N+2]=log(c);
        
        
      }
    }
    pick.th=rbinom(1,1,.5);
    
    temp[2*N+7]=z[2*N+7]+zn[2*N+3];
    
    u=log(runif(1,0,1));
    alph=logf1(2,temp,E1.old)-logf1(2,z,E1.old)-temp[2*N+7]+z[2*N+7];
    if (alph>=u) {
      z[2*N+7]=temp[2*N+7];
      acc[5]=acc[5]+1;
    }
    
    
    
    if(i %% 1000==0) {
      if(i <=burn ) {
        c0[1:5]=c0[1:5]+(acc[1:5]/1000>.2)*.5*c0[1:5]-(acc[1:5]/1000<.15)*.5*c0[1:5];
        c0[6:7]=c0[6:7]-(acc[6:7]/1000>.4)*c0[6:7]*(1-c0[6:7])+(acc[6:7]/1000<.1)*c0[6:7]*(1-c0[6:7]);
        c0[6:7]=sapply(c0[6:7],function(x) max(.6,min(.99,x)));
      }
      print(acc/1000);
      acc=rep(0,7);	
    }
    MC.tot[i,]=z;
  }
  print(l);
  MC=MC.tot[-(1:burn),];
  d=density(MC[,2*N+1]);
  if (l==1) {
    plot(d,xlim=c(0,1),main="Posterior Density");
  } else {
    lines(d,col=l);
  }
  plot(density(MC[,2*N+4]),xlim=c(0,1),main="Posterior Density");
  rho1[,l]=MC[,2*N+1]; rho2[,l]=MC[,2*N+4];
  G1[,l]=MC[,2*N+3]; G2[,1]=MC[,2*N+6];
  
  for (v in 1:10) {
    d1=density(MC[,2*N+1]);
    max.d1[l]=d$x[which.max(d1$y)];
    d2=density(MC[,2*N+4]);
    max.d2[l]=d$x[which.max(d2$y)];
  }
  P1[l]=mean(post.pi(MC[,2*N+1]));
  P2[l]=mean(post.pi(MC[,2*N+4]));
  

  
  
}
stop_time=Sys.time()
stop_time-start_time;

post.pi=function(rho) {
  return(dbeta(rho,alpha,beta)/(dbeta(rho,alpha,beta)+dbeta(rho,alpha1,beta1)));
}
post.pi=Vectorize(post.pi)

