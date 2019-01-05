install.packages("gsl")
install.packages("PearsonDS")
library(gsl)
library(PearsonDS)
?PearsonDS
qpearson(0.5,moments = c(10.5,3.142^2,1.1,5.6))
qpearson(0.99865,moments = c(10.5,3.142^2,1.1,5.6))
qpearson(0.00135,moments = c(10.5,3.142^2,1.1,2.6))
test=pearsonFitML(10.5,3.142^2,1.1,2.6)
qpearsonI(p=0.5,a=0.1346417,b=0.393575,location=8.228176,scale=8.912659)

aranesp=read.csv("normal_stats_output_all_Oct5.csv")
View(aranesp)
