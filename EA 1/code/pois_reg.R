### Poisson Regression ####
source("code/library.R")
source("code/data_format.R")
source("code/helpers.R")

# fitted poisson, loaded for faster paper running
fit.pois <- dlply(dat, .(store), function(x) {
  glm(data=x, mvm ~ price, family=poisson, maxit=50)
})

save(fit.pois, file="fit.pois.rda")

