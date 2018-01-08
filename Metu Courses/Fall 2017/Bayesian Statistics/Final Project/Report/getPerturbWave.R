getPerturbWave<-function(A,A_0,x_0,x,sig_inc,sig_top,N)
{

A_inc<-A_0*exp(-0.5*x^2/sig_inc^2) #incident wave field
h<-A*exp(-0.5*(x)^2/sig_top^2)
y<-log(A_inc)+log(h)+rnorm(N,0,5)

return(y)
}