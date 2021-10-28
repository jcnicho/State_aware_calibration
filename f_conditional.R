logf1=function(k,z,E) {
  z[2*N+2]=exp(z[2*N+2]); z[2*N+5]=exp(z[2*N+5]); z[2*N+7]=exp(z[2*N+7])
  if (k==1) { #rho1
    bz=diag(1/sqrt(abs(E$values)))%*%t(E$vectors)%*%as.matrix(z[1:N]);
    result=-1/2*sum(log(abs(E$values)))-(e/2+N/2)*log(1+sum(bz^2));
    result=result+z[2*N+3]*dbeta(z[2*N+1],alpha,beta,log=TRUE)+(1-z[2*N+3])*dbeta(z[2*N+1],alpha1,beta1,log=TRUE);
  } else if (k==2) { #sig2
    result=sum(dnorm(eta(x,z[1:N],z[(N+1):(2*N)]),mean=y,sd=sqrt(z[2*N+7]),log=TRUE));
    result=result+dinvgamma(z[2*N+7],shape=1/w,scale=w,log=TRUE);
    
  } else if (k==4) { #theta1
    result=sum(dnorm(eta(x,z[1:N],z[(N+1):(2*N)]),mean=y,sd=sqrt(z[2*N+7]),log=TRUE));
    bz=diag(1/sqrt(abs(E$values)))%*%t(E$vectors)%*%as.matrix(z[1:N]);
    result=result-1/2*sum(log(abs(E$values)))-(e/2+N/2)*log(1+sum(bz^2));
  } else if (k==5) { #rho2
    bz=diag(1/sqrt(abs(E$values)))%*%t(E$vectors)%*%as.matrix(z[(N+1):(2*N)]);
    result=-1/2*sum(log(abs(E$values)))-(e/2+N/2)*log(1+sum(bz^2));
    result=result+z[2*N+6]*dbeta(z[2*N+4],alpha,beta,log=TRUE)+(1-z[2*N+6])*dbeta(z[2*N+4],alpha1,beta1,log=TRUE);
  } else { #theta2
    result=sum(dnorm(eta(x,z[1:N],z[(N+1):(2*N)]),mean=y,sd=sqrt(z[2*N+7]),log=TRUE));
    bz=diag(1/sqrt(abs(E$values)))%*%t(E$vectors)%*%as.matrix(z[(N+1):(2*N)]);
    result=result-1/2*sum(log(abs(E$values)))-(e/2+N/2)*log(1+sum(bz^2));
  }
  return(result);
}