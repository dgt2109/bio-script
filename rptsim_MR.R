
#params = c(2,4,6,8,10)
#U: params = c(0.5,5,10,20)
#oY: params = c(10,20,30,40)

parameters = c(10,20,30,40,50)
allbetas = list()
allpower = matrix(0,length(parameters),5)

for (p in 1:length(parameters)) {
  param = parameters[p]
  source("simulate_MR.R")
  allpower[p,] = power
  allbetas[[p]] = betas
}

write.csv(allpower,"C:/Users/david/Desktop/tp_power.csv")
write.csv(allbetas,"C:/Users/david/Desktop/tp_betas.csv")

