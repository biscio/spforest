expect_silent(arbre <- intensitytree(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  mtry = 1,
  minpts = 500
))

expect_equal(names(arbre), 
             c("tree", "X", "namecov", 
               "namelist", "listcov", "im"))

expect_length(arbre$namecov, 15)


expect_length(arbre$tree, 23)

expect_true(spatstat.geom::is.im(arbre$im))

expect_true(spatstat.geom::is.ppp(arbre$X))


expect_equal(arbre$tree[[4]]$nodeID, 4)

expect_equal(arbre$tree[[4]]$left_daughter, 8) 
expect_equal(arbre$tree[[4]]$right_daughter, 9)  

expect_equal(arbre$tree[[4]]$split_var, 15)

expect_equal(arbre$tree[[4]]$split_val, 4.599173,
             tolerance = 1e-5)

expect_equal(arbre$tree[[4]]$status, 1)

expect_equal(arbre$tree[[4]]$intensity_pred, 0.00862606,
             tolerance = 1e-5
)


A <- intensitytree(
  X = spatstat.data::bei,
  listcovariates = beisoilres,
  mtry = 1,
  minpts = 100
)
B <- format(object.size(A), units = "Mb")
expect_true(as.numeric(gsub(" Mb", "",B)) < 4)

# 
# ftree <- function() {
#   treerec(
#     X = spatstat.data::bei,
#     listcovariates = beisoilres,
#     mtry = 2 / 3,
#     minpts = 20
#   )
# }
# 
# gtree <- function() {
#   intensitytree(
#     X = spatstat.data::bei,
#     listcovariates = beisoilres,
#     mtry = 1,
#     minpts = 100
#   )
# }
# 
# 
# 
# library(microbenchmark)
# microbenchmark(ftree(), gtree(), times = 10)
