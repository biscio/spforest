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


# Generate thicker resolution for vignettes and testings purposes ----

beisoilsmall <- lapply(beisoilres, FUN=function(i){
  as.im(i, dimyx=c(10,20))
})

plot(beisoilsmall[[8]])
