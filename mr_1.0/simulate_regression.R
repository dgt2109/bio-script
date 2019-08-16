# model: y = m (d + n) + e

niter = 300
nsucc = 0

for (i in 1:niter) {
    
  d1 = 1
  d2 = 5
  m = 0
  
  eta1 = runif(400,-1,1)
  eta2= runif(5,-1,1)
  
  eps1 = rnorm(length(eta1),0,2)
  eps2 = rnorm(length(eta2),0,2)
  
  x = c(d1+eta1, d2+eta2)
  eps = c(eps1,eps2)
  y = m * x + eps
  
  plot(x,y)
  res = summary(lm(y~x))
  if (res$coefficients[8] <0.05) {
    nsucc = nsucc + 1
  }
}
print(nsucc/niter)