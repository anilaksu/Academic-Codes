####################################################
#	                                             #
#      Posterior, Perspective and Contour plots    #
#                  by Anil Aksu                    #    
#     						         #
####################################################	

## the range of sampling
x=seq(0,20,length=101)
## this function gets numbers from console
posterior=dnorm(x, mean = 7, sd = 1.5, log = FALSE)


## let's plot them
plot(range(x), range(c(posterior)), type='n', xlab=expression(paste(theta)), ylab=expression(paste("f(", theta, "|x )")))

lines(x, posterior, type='l', col='blue',lwd=5)

title("Posterior Distribution")
  legend = c("posterior")

## perspective plot
x <- seq(-10, 10, length= 30)
y <- x
f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, f)
z[is.na(z)] <- 1
op <- par(bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue")
persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "lightblue",
      ltheta = 120, shade = 0.75, ticktype = "detailed",
      xlab = "X", ylab = "Y", zlab = "Sinc( r )"
) -> res
round(res, 3)

# contour plot
a <- expand.grid(1:20, 1:20)
b <- matrix(a[,1]^2 + a[,2]^2, 20)
filled.contour(x = 1:20, y = 1:20, z = b,
               plot.axes = { axis(1); axis(2); points(10, 10) })


## bivariate posterior sampling

## the range of sampling
x=seq(-4,6,length=101)
## this function gets numbers from console
posterior=0.8*dnorm(x, mean = 0, sd = 1, log = FALSE)+0.2*dnorm(x, mean = 4, sd = 1, log = FALSE)

## let's plot them
plot(range(x), range(c(posterior)), type='n', xlab=expression(paste(theta)), ylab=expression(paste("f(", theta, "|x )")))

lines(x, posterior, type='l', col='blue',lwd=5)

# title("Bivariate Posterior Distribution")
  legend = c("posterior")


