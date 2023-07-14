library(profvis)

Xtest

target.points

profvis(tesstree(X=Xtest, lambda=100, 
                 target.points = target.points, test.connected = T))



a<-tessforest(X=Xtest, lambda=NULL, Ntree=100, mc.cores = 5, 
            at = NULL, test.connected = T)
plot(a)
