#-------------------Function Script: Pearson Ppk Calculation Pipeline
#Input: df containing 4 moments and other aggregate stats for each parameter
#Output: same as input, with a "Pearson.Ppk" column added; output file name: input's file name
#appended with "_pearson"
#Created on: Oct 22


setwd("~/OliviaDai")
install.packages("gsl")
install.packages("PearsonDS")


library(gsl)
library(PearsonDS)
print(noquote("Pearson Ppk Calculation Program"))
print(noquote("Put your dataset under the working directory."))
df_path=readline(prompt = "Enter file name(Don't include extention):")
df_path_no_OOS=readline(prompt = "Enter no_OOS file name(Don't include extention):")
#eg:df_path=program+"summary_stats_pearson.csv"
#           program+"summary_stats_no_OOS_pearson.csv"

Pearson_ppk=function(df, input_type){
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
    if (is.na(df[i,"spec status"]))
    {
      ppk=NA
    }
    if(!is.na(df[i,"spec status"]))
    {
      if (df[i,"spec status"]==0){
        ppk=min(c((UAL-m50)/(Up-m50),(m50-LAL)/(m50-Lp)))
      }
      else if (df[i,"spec status"]==1){
        ppk=(UAL-m50)/(Up-m50)
      }
      else if (df[i,"spec status"]==-1){
        ppk=(m50-LAL)/(m50-Lp)
      }
      
    }
    if (input_type=="all")
      {
        df[i,"pearson Ppk"]=ppk
    }
    else if(input_type=="no_OOS")
    {
      df[i,"pearson Ppk_no_OOS"]=ppk
    }
    
    
  }
  return(df)
  
}

#input_df: the df containing Aggre. stats for each parameter;used for pearson ppk calculation;
#columns it must have:
#"parameter name",'unit procedure',(N),"mean","median","STD","variance","UAL","LAL",
#"spec status",(Ppk_noOOST),"skewness","kurtosis"
#saving functions as a workspace or Rdata file?
#input_type: "all" OR "no_OOS" OR "no_OOST"...[TO BE ADDRESSED]

Pearson_data_handling=function(df_path,input_type)
{
  input_df=read.csv(paste(df_path,".csv",sep=""),check.names=FALSE)
  #create a new column to store the new output
  if (input_type=="all")
    {
    input_df["pearson Ppk"]=NA
    }
  else if (input_type=="no_OOS")
    {
      input_df["pearson Ppk_no_OOS"]=NA
    }
  #handle the spec format issue: set missing spec to NA
  input_df[input_df["UAL"]==".","UAL"]=NA
  input_df[input_df["LAL"]==".","LAL"]=NA
  
  #transform the string-typed specs to Numerical
  input_df["UAL"]=lapply(lapply(input_df["UAL"],as.character),as.numeric)
  input_df["LAL"]=lapply(lapply(input_df["LAL"],as.character),as.numeric)
  
  #---------Dataset preparation complished, send the df to Pearson calculation func
  output_df=Pearson_ppk(input_df, input_type)
  
  write.csv(output_df,file = paste(df_path,"pearson.csv",sep="_"),row.names = FALSE)
  
}

# Implement the function
Pearson_data_handling(df_path, input_type = "")

Pearson_data_handling(df_path_no_OOS, input_type = "no_OOS")

#test eg: test=read.csv("PearsonPpk_normal_stats_output_all_noOOST_Oct22.csv")
