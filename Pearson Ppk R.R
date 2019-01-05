install.packages("gsl")
install.packages("PearsonDS")
setwd("~/OliviaDai")
library(gsl)
library(PearsonDS)
?PearsonDS
qpearson(0.5,moments = c(10.5,3.142^2,1.1,5.6))
qpearson(0.99865,moments = c(10.5,3.142^2,1.1,5.6))
qpearson(0.00135,moments = c(10.5,3.142^2,1.1,2.6))
test=pearsonFitML(10.5,3.142^2,1.1,2.6)
qpearsonI(p=0.5,a=0.1346417,b=0.393575,location=8.228176,scale=8.912659)


Pearson_ppk=function(df){
  N=dim(df)[1]
  
  for (i in 1:N){
    #kurtosis adjustment
    single_row=df[i,]
    moments_input=df[i,c("mean","variance","skewness","kurtosis")]+c(0,0,0,3)
    
    m50=qpearson(0.5,moments = moments_input)
    Lp=qpearson(0.00135,moments = moments_input)
    Up=qpearson(0.99865,moments = moments_input)
    UAL=single_row[["UAL"]]
    LAL=single_row[["LAL"]]
    if (is.na(df[i,"spec.status"]))
    {
      ppk=NA
    }
    if(!is.na(df[i,"spec.status"]))
      {
      if (df[i,"spec.status"]==0){
        ppk=min(c((UAL-m50)/(Up-m50),(m50-LAL)/(m50-Lp)))
      }
      else if (df[i,"spec.status"]==1){
        ppk=(UAL-m50)/(Up-m50)
      }
      else if (df[i,"spec.status"]==-1){
        ppk=(m50-LAL)/(m50-Lp)
      }
    
    }
    df[i,"pearson.ppk"]=ppk

  }
  return(df)
  
}
#Note: the SPACE in the column names are auto replaced by . in R, UAL and LAL are read as factors
#R functions do not change the original copies of dfs, need to return a new df
aranesp=read.csv("normal_stats_output_all_Oct11.csv")
aranesp_no_oost=read.csv("normal_stats_output_all_noOOST_Oct22.csv")

aranesp["pearson.ppk"]=NA
aranesp_no_oost["pearson.ppk"]=NA

View(aranesp)
View(aranesp_no_oost)

aranesp[aranesp["UAL"]==".","UAL"]=NA
aranesp[aranesp["LAL"]==".","LAL"]=NA

aranesp_no_oost[aranesp_no_oost["UAL"]==".","UAL"]=NA
aranesp_no_oost[aranesp_no_oost["LAL"]==".","LAL"]=NA

aranesp["UAL"]=lapply(lapply(aranesp["UAL"],as.character),as.numeric)
aranesp["LAL"]=lapply(lapply(aranesp["LAL"],as.character),as.numeric)

aranesp_no_oost["UAL"]=lapply(lapply(aranesp_no_oost["UAL"],as.character),as.numeric)
aranesp_no_oost["LAL"]=lapply(lapply(aranesp_no_oost["LAL"],as.character),as.numeric)

aranesp_pearsonppk=Pearson_ppk(aranesp)
aranesp_no_oost_pearsonppk=Pearson_ppk(aranesp_no_oost)


View(aranesp_pearsonppk)
View(aranesp_no_oost_pearsonppk)
write.csv(aranesp_pearsonppk,file = "aranesp_pearsonppk_stats_Oct11.csv")
write.csv(aranesp_no_oost_pearsonppk,file = "aranesp_no_oost_pearsonppk_stats_Oct22.csv",
          row.names = FALSE)

