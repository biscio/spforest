# I have forgotten where I got the soils data.
# beisoil <- list(
#     bei.extra$elev, bei.extra$grad,
#     Al, B, Ca, Cu, Fe, K, Mg, Mn, P, Zn, N, Nmin, pH
#   )

# Putting all the data to the same resolution as bei.extra$elev ----

beisoilres <- NULL
beisoilres[[1]] <- beisoil[[1]]
beisoilres[[2]] <- beisoil[[2]]
for (i in 3:15) {
  W <- commonGrid(beisoil[[1]], beisoil[[i]])
  beisoilres[[i]] <- beisoil[[i]][W, drop = F]
  rm(W)
}
names(beisoilres) <- c(
  "elev", "grad", "Al", "B", "Ca",
  "Cu", "Fe", "K", "Mg", "Mn",
  "P", "Zn", "N", "Nmin", "pH"
)

# Scale and normalised all images ----

beisoilnorm <- lapply(beisoilres, FUN = function(i) {
  (i - mean(i)) / sd(i)
})

names(beisoilnorm) <- c(
  "elev", "grad", "Al", "B", "Ca",
  "Cu", "Fe", "K", "Mg", "Mn",
  "P", "Zn", "N", "Nmin", "pH"
)


# Downscaled version of the data for vignettes and testing purposes ----

beisoilxsmall <- lapply(beisoilres, FUN=function(i){
  as.im(i, dimyx=c(10,20))
})

beisoilxsmall <- lapply(beisoilres, FUN=function(i){
  as.im(i, dimyx=c(5,10))
})

usethis::use_data(beisoilsmall, compress = "xz")
usethis::use_data(beisoilxsmall, compress = "xz")

# save(beisoilsmall, 
#      file = "data/beisoilsmall.RData", 
#      compress = "xz")
# 
# save(beisoilxsmall, 
#      file = "data/beisoilxsmall.RData", 
#      compress = "xz")

plot(beisoilsmall[[8]])
