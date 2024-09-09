# TODO: redo with new version splitcell

# 
# expect_silent(
#   output <- splitcell(
#     X = spatstat.data::bei[seq.int(from = 1, to = 3604, by = 100)],
#     score = "lcv",
#     listcovariates = list(
#       grad = spatstat.data::bei.extra$grad,
#       elev = spatstat.data::bei.extra$elev
#     ),
#     usecovariates = c(1, 1),
#     thres.cell = 100,
#     minpts = 10,
#     tol = Inf,
#     imp = NULL,
#     minsplitq = 0.5,
#     maxsplitq = 0.5
#   )
# )
# 
# expect_equal(spatstat.geom::npoints(output$PPleft), 14)
# 
# expect_equal(spatstat.geom::npoints(output$PPright), 23)
# 
# expect_equal(output$split_val, 0.06202988)
# 
