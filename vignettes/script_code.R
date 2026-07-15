library(mahiMSE)

Hist = Simulate(smallOM)

StatQuo = HalfEff = mahiMP
formals(HalfEff)$rel_TAE = 0.5
class(StatQuo) = class(HalfEff) = 'mp'

myMSE = Project(Hist, c('StatQuo','HalfEff'))

mahiplot(myMSE)
