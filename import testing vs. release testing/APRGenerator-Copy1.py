import pandas as pd
import numpy as np

import os
import re

import matplotlib.pyplot as plt  

from scipy.stats import anderson
from scipy.stats import kstest
from scipy.stats import shapiro
#Nov 6: examined and worked well. Further update the spec correction function to groupby
#no need to do the conversion if the df is read from excel. Do for csv
#two string_to_num functions: both require string-type input
#update on Nov: these two conversion functions do not address point in float. MARKED DOWN
#Notes on Dec 6:
#For Enbrel, the specs may only contain text contents,specs convertion function will cause error.
#Solution: consider first limiting to quantitative data and then converting the specs.

#Dec 7: bug found in spec_correction_rawset---this function needs to be modified if the grouping changes(eg: by prds, mtd, sites..)
        #specs may be populated in a wrong way if using the wrong grouping(multiple methods per combination)
#Dec 10: consider ctrl_num_conversion as well

#Dec 18: modified para_unit_combi_df--allow for flexible grouping as input. return unique combinations of columns specified.

def string_to_num(input_cell):
    
    
    #update on Nov 14: add judgement on input type. directly output the numerial values. 
    #now able to handle a mixture of int/float and str
    if type(input_cell)!=str:
        #compare to str because the specs can be float or int
        return(input_cell)
    
    if (input_cell=="."):
        return(".")
    else:
        #address the non-missing and string-typed numerical specs
        newstr = ''.join((ch if ch in '0123456789.' else ' ') for ch in str(input_cell))
        result_list=[float(i) for i in newstr.split()]
        #Dec 7:
        if (result_list==[]):
            #pure textual spec
            return(input_cell)
        else:
            
            return(result_list[0])


        
class APR_analytical_pipeline:
    def __init__(self, program, file_path ):
        self.program=program
        self.df=pd.read_csv(file_path,index_col=False)
        
        
    def kurt(self, values):
        n=len(values)
        if (n<4): 
            result=np.nan
        else:
            sample_mean=np.mean(values)
            sample_std=np.std(values,ddof=1)
            result=np.sum(((values-sample_mean)/sample_std)**4)*n*(n+1)/((n-1)*(n-2)*(n-3))-3*(n-1)**2/((n-2)*(n-3))
        return(result)

    def skewness(self, values):
        n=len(values)
        if (n<3):
            result=np.nan
        else:
            sample_mean=np.mean(values)
            sample_std=np.std(values,ddof=1)
            result=n*np.sum(((values-sample_mean)/sample_std)**3)/((n-1)*(n-2))
        return(result)

    def CV(self, values):
        sample_mean=np.mean(values)
        sample_std=np.std(values,ddof=1)
        return(100*sample_std/sample_mean)
    
    def shapiro_func(self,values):
        if (len(values)<4):
            return(np.nan)
        else:
            return(shapiro(values)[1])
              
         
    def kstest_func(self, values):
        if len(values<4):
            return(np.nan)
        else:
            return(kstest(values, "norm")[1])

    def ppk_single_row(self, values,ual,lal,spec_status):
        if spec_status==0:
            ppk=np.minimum(float(ual)-np.mean(values),np.mean(values)-float(lal) )/(3*np.std(values,ddof=1))
        else:
            if spec_status==1:
                ppk=(float(ual)-np.mean(values))/(3*np.std(values,ddof=1))
            else:
                if spec_status==-1:
                    ppk=(np.mean(values)-float(lal))/(3*np.std(values,ddof=1))
                else:
                    ppk="NA3"
        return(ppk)

    def Ppk_groupby(self, df, min_size):
        #try change the index to -1 so as to avoid mistakes from specification correction 
        spec_status=df["spec status"].iloc[0]
        ual=df.UAL.iloc[0]
        lal=df.LAL.iloc[0]
        values=df['value_num(recorded/full precision)']

        if (len(values)<min_size):
                    ppk="NA1"

        if( (len(values)>=min_size) and (len(pd.unique(values)<5)) ):
                    ppk="NA2"

        if( (len(values)>=min_size) and (len(pd.unique(values)>=5)) ):

            if spec_status==0:
                ppk=np.minimum(float(ual)-np.mean(values),np.mean(values)-float(lal) )/(3*np.std(values,ddof=1))
            else:
                if spec_status==1:
                    ppk=(float(ual)-np.mean(values))/(3*np.std(values,ddof=1))
                else:
                    if spec_status==-1:
                        ppk=(np.mean(values)-float(lal))/(3*np.std(values,ddof=1))
                    else:
                        ppk="NA3"
        #return a single value will make it similar to std, mean function, easy for group by but hard for generating summary tables
        #can also do: df[ppk]=ppk, return(df), the df will be changed outside the function. same ppk populated for each group
        #Raw ppk(with oos)
        df["Ppk"]=ppk
        return(ppk)
    
     #Nov 14: bug resolved for comparing a numerical series to "."     
     #Indices need to be reset before using these functions because the map function disturb the original setup
    #no need to reset if we use df.UAL==... but this code does not work for pure numerical series
    def add_specside(self,df):
        #no spec default
        df["spec status"]=None
        
        #two specs
        df["spec status"].loc[(pd.Series(map(str,df.UAL))!=".")&(pd.Series(map(str,df.LAL))!=".")]=0
        #df["spec status"].loc[(df.UAL!=".") & (df.LAL!=".")]=0
        
        #only lower spec
        df["spec status"].loc[(pd.Series(map(str,df.UAL))==".")&(pd.Series(map(str,df.LAL))!=".")]=-1
        #df["spec status"].loc[(df.UAL==".") & (df.LAL!=".")]=-1
        
        #only upper spec
        df["spec status"].loc[(pd.Series(map(str,df.UAL))!=".")&(pd.Series(map(str,df.LAL))==".")]=1
        #df["spec status"].loc[(df.UAL!=".") & (df.LAL==".")]=1
        
        return(df)

    def add_clside(self, df):
        #no spec default
        df["ctrl status"]=None
        
        #two specs
        df["ctrl status"].loc[(pd.Series(map(str,df.UCL))!=".")&(pd.Series(map(str,df.LCL))!=".")]=0
        #df["ctrl status"].loc[(df.LCL!=".") & (df.UCL!=".")]=0
        
        #only lower spec
        df["ctrl status"].loc[(pd.Series(map(str,df.UCL))==".")&(pd.Series(map(str,df.LCL))!=".")]=-1
        #df["ctrl status"].loc[(df.UCL==".") & (df.LCL!=".")]=-1
        #only upper spec
        df["ctrl status"].loc[(pd.Series(map(str,df.UCL))!=".")&(pd.Series(map(str,df.LCL))==".")]=1
        #df["ctrl status"].loc[(df.UCL!=".") & (df.LCL==".")]=1

        return(df)
    
    def df_colname_unify(self,df):
        
        df.columns = map(str.lower, df.columns)
        df=df.rename(columns={"ual":"UAL","lal":"LAL","lcl":"LCL","ucl":"UCL","result_type":"result type",\
                             'value_num(recorded full precision)':'value_num(recorded/full precision)'})
        
        return(df)
    
    def specs_num_conversion(self,df):
        
        #Excel file is able to distinguish string and int in the same column. 
        #No need to run this func for a df read-in from excel(will cause error inside the sting_to_num function.
        
        #Do run this if the df is read from a csv file

        
        #include a return so as to change the df outside the function space, write df=function() to update in the main
        df["UAL"]=pd.Series(map(string_to_num,df["UAL"]))
        df["LAL"]=pd.Series(map(string_to_num,df["LAL"]))
        
        return(df)
        
    def quanti_quali_namelists(self,df):
        
        
        qualitative_paras=list(pd.unique(df["parameter name"][pd.Series(map(str.lower,df['result type']))=="text"]))

        quantitative_paras=list(pd.unique(df["parameter name"][pd.Series(map(str.lower,df['result type']))=="number"]))

        return(quantitative_paras,qualitative_paras)

    #Maintain flexible input, allow for df, df_quan and other dataframe input 
    def para_unit_combi_df(self, df, grouping_list=["parameter name","unit procedure"]):

        para_unit_combi=df[df.columns.intersection(grouping_list)].copy().drop_duplicates().reset_index(drop=True)
        para_unit_combi=para_unit_combi.sort_values(by=grouping_list)
        #only return result, do not create new variable in the class
        return(para_unit_combi)

    #only use a quantitative dataset as the input
    #NOTE: this function needs to be modified if the grouping changes(eg: by prds, mtd, sites..)
    def spec_correction_rawset(self, df_quan):
        
        para_unit_pair_df=self.para_unit_combi_df(df_quan)
        
        for i in range(len(para_unit_pair_df)):
            para=para_unit_pair_df["parameter name"].iloc[i]
            unit=para_unit_pair_df["unit procedure"].iloc[i]

            sub_df_len=len(df_quan.loc[(df_quan['parameter name']==para) &(df_quan["unit procedure"]==unit)])

            #more than one row, needs spec correction
            if (sub_df_len>1):

                #could loop using the uniq combination of unit and para, otherwise add this row number check to avoid void error
                ual=df_quan[(df_quan['parameter name']==para) &(df_quan["unit procedure"]==unit)].UAL.iloc[-1]

                lal=df_quan[(df_quan['parameter name']==para) &(df_quan["unit procedure"]==unit)].LAL.iloc[-1]
                df_quan["UAL"].loc[(df_quan['parameter name']==para) &(df_quan["unit procedure"]==unit)]=ual
                df_quan["LAL"].loc[(df_quan['parameter name']==para) &(df_quan["unit procedure"]==unit)]=lal

        return(df_quan)   
     
        #ok to remove coerce to numeric if we run the spec_number_conversion first;
        #only need this coersion if the df is directly read from csv
    def find_oos_index(self, df):
        df_twospecs=df[df["spec status"]==0]
        df_uspecs=df[df["spec status"]==1]
        df_lspecs=df[df["spec status"]==-1]

        #Find the indexes of OOS rows, try to do this in SQL
        twospecs_rows_oos=list(df_twospecs.index[(df_twospecs["value_num(recorded/full precision)"]<\
                                                  pd.to_numeric(df_twospecs.LAL, errors="coerce"))|\
                        (df_twospecs["value_num(recorded/full precision)"]>pd.to_numeric(df_twospecs.UAL, errors="coerce"))])

        uspecs_rows_oos=list(df_uspecs.index[(df_uspecs["value_num(recorded/full precision)"]> \
                                              pd.to_numeric(df_uspecs.UAL,errors="coerce"))])

        lspecs_rows_oos=list(df_lspecs.index[(df_lspecs["value_num(recorded/full precision)"]<\
                                              pd.to_numeric(df_lspecs.LAL,errors="coerce"))])

        return(list(set(twospecs_rows_oos+uspecs_rows_oos+lspecs_rows_oos)))

    def find_oot_index(self, df):
        df_twocls=df[df["ctrl status"]==0]
        df_ucl=df[df["ctrl status"]==1]
        df_lcl=df[df["ctrl status"]==-1]

        twocls_rows_oot=list(df_twocls.index[(df_twocls["value_num(recorded/full precision)"]<\
                                              pd.to_numeric(df_twocls.LCL,errors="coerce"))\
                    |(df_twocls["value_num(recorded/full precision)"]>pd.to_numeric(df_twocls.UCL, errors="coerce"))])

        ucls_rows_oot=list(df_ucl.index[(df_ucl["value_num(recorded/full precision)"]>pd.to_numeric(df_ucl.UCL,errors="coerce"))])

        lcls_rows_oot=list(df_lcl.index[(df_lcl["value_num(recorded/full precision)"]<pd.to_numeric(df_lcl.LCL,errors="coerce"))])

        return(list(set(twocls_rows_oot+ucls_rows_oot+lcls_rows_oot)))
                   
                   
    def aggre_summary_df_generator(self, df, input_type, min_size):
        
        #input_type:"" OR "_no_OOS"
        
        #initialize the summary output df
        #aggre_stats_output_df=pd.DataFrame(columns=["parameter name",'unit procedure',#"PRDS",#\
                                              # "N","mean","median","STD","variance","UAL","LAL","spec status","Ppk"+input_type,\
                                             #  "skewness","kurtosis","CV",#"anderson-darling"\
                                               # "shapiro","KS"\
                                                #  ])
        
        aggre_stats_output_df=pd.DataFrame(data={\
                "N":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].count(),\
              "mean":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].mean(),\
          "median":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].median(),\
          "STD":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].std(),\
          "varaince":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].var(),\
          "Ppk"+input_type:df.groupby(["parameter name","unit procedure"]).apply(self.Ppk_groupby,min_size=min_size),\
          "skewness":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].apply(self.skewness),\
          "kurtosis":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].apply(self.kurt),\
          "CV":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].apply(self.CV),\
          "shapiro":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].apply(self.shapiro_func),\
          "KS":df.groupby(["parameter name","unit procedure"])["value_num(recorded/full precision)"].apply(self.kstest_func),\
          "spec status":df.groupby(["parameter name","unit procedure"])["spec status"].last(),\
          "UAL":df.groupby(["parameter name","unit procedure"])["UAL"].last(),\
          "LAL":df.groupby(["parameter name","unit procedure"])["LAL"].last()\
                                            })
        return(aggre_stats_output_df)
        
        
        
        
        


    #change the threshold(30) if we allow smaller sample size after excluding the OOS
    #input_type can be "" or "_no_OOS" or "no_OOST", etc, used to define the the column name
    
    def aggre_summary_df_generator_old(self, df, input_type, min_size=30):
        
        #input_type:"" OR "_no_OOS"
        
        #initialize the summary output df
        aggre_stats_output_df=pd.DataFrame(columns=["parameter name",'unit procedure',#"PRDS",#\
                                               "N","mean","median","STD","variance","UAL","LAL","spec status","Ppk"+input_type,\
                                               "skewness","kurtosis","CV",#"anderson-darling"\
                                                "shapiro","KS"\
                                                  ])
        
        para_unit_pair_df=self.para_unit_combi_df(df)
        
        for i in range(len(para_unit_pair_df)):

            para=para_unit_pair_df["parameter name"].iloc[i]
            unit=para_unit_pair_df["unit procedure"].iloc[i]

            single_df=df.loc[(df['parameter name']==para) &(df["unit procedure"]==unit)].copy().sort_values(by="parameter date") 
            single_df_values=single_df["value_num(recorded/full precision)"]  

            ual=single_df.UAL.iloc[-1]
            lal=single_df.LAL.iloc[-1]
            spec_status=single_df["spec status"].iloc[-1]

            if (len(single_df)<min_size):
                ppk="NA1"

            if( (len(single_df)>=min_size) and (len(pd.unique(single_df_values)<5)) ):
                ppk="NA2"

            #NA3: no specs; addressed inside the ppk_single_row

            if ((len(single_df)>=min_size) and (len(pd.unique(single_df_values)>=5)) ):

                #prds=pd.unique(single_df["parameter detail"])

                ppk=self.ppk_single_row(single_df_values,ual,lal,spec_status)
                
                
            if (len(single_df)<4):
                shap=np.nan
                ks=np.nan
            else:
                shap=shapiro(single_df_values)[1]
                ks=kstest(single_df_values,"norm")[1]

            aggre_stats_output_df=aggre_stats_output_df.append({\
                            'parameter name':para,\
                            "unit procedure":unit,\
                            #"PRDS":prds,
                            "N":len(single_df_values),\
                            "mean":np.mean(single_df_values),\
                            "STD":np.std(single_df_values,ddof=1),\
                            "variance":np.var(single_df_values,ddof=1),\
                            "median":np.median(single_df_values),\
                            "LAL":lal, "UAL":ual,      \
                            "spec status":spec_status ,\
                            "Ppk"+input_type: ppk,\
                            "skewness":self.skewness(single_df_values),\
                            "kurtosis":self.kurt(single_df_values),\
                            "CV":self.CV(single_df_values),\
                            #"Anderson-Darling":anderson(sub_df),\
                            "shapiro":shap,\
                            #the p value of shapiro
                           "KS":ks},\
                            #the p value of kstest
                                        ignore_index=True)
        return(aggre_stats_output_df)
    
    def summary_print(self, df, purpose,input_type):
        #R_pearson: remove all rows with a NA on the ppk
        #Report: just output the entire df without changing anything(keep NAs on Ppk)
        #input_type: "" OR "_no_OOS"
        if (purpose=="R_pearson"):
            if (input_type=="_no_OOS"):
                df[(df["Ppk_no_OOS"]!="NA1")&(df["Ppk_no_OOS"]!="NA2")&(df["Ppk_no_OOS"]!="NA3")].\
                to_csv(self.program+"_summary_stats"+input_type+".csv",index=False)
                
            if (input_type==""):
                df[(df["Ppk"]!="NA1")&(df["Ppk"]!="NA2")&(df["Ppk"]!="NA3")].to_csv(self.program+"_summary_stats"+input_type+".csv",index=False)

        
        if(purpose=="Report"): 
            df.to_csv(self.program+"_summary_stats"+input_type+"_Report"+".csv", index=False)
          
        return

    def JMP_single_para_subset_generator(self, extract_df , raw_df, folder_path):
        
        para_unit_pair_df=self.para_unit_combi_df(extract_df)
        os.chdir(folder_path)
        
        for i in range(len(para_unit_pair_df)):
            para=para_unit_pair_df["parameter name"].iloc[i]
            unit=para_unit_pair_df["unit procedure"].iloc[i]
            single_subset=raw_df[(raw_df["parameter name"]==para)&(raw_df["unit procedure"]==unit)].copy()
            single_subset.to_csv(para+"-"+unit+"JMP.csv", index=False)

        return
    
    def qqplot_screening_shapiro(self, program,shapiro_alpha,raw_df,summary_df):
        raw_summary_merged=raw_df.merge(summary_df,on=["parameter name","unit procedure"])
        para_unit_comb=raw_summary_merged[(raw_summary_merged["shapiro"]<shapiro_alpha)&\
                                          (~raw_summary_merged["Ppk_no_OOS"].isin(["NA1","NA2","NA3"]))&(raw_summary_merged["N"]>=30)]\
                                                        [["parameter name","unit procedure"]].drop_duplicates()

        fig=plt.figure(figsize=(20,24))
        fig.suptitle(program+"\n"+"Q-Q plot:"+"Shapiro<"+str(shapiro_alpha),fontsize=20)

        for i in range(len(para_unit_comb)):
            ax=plt.subplot(int(len(para_unit_comb)/5)+1,5,i+1,frameon=False)

            stats.probplot(raw_summary_merged["value_num(recorded/full precision)"][(raw_summary_merged["parameter name"]==\
                             para_unit_comb["parameter name"].iloc[i])&\
                            (raw_summary_merged["unit procedure"]==para_unit_comb["unit procedure"].iloc[i])],\
                          plot=plt);
            ax.set_title(str(para_unit_comb["parameter name"].iloc[i])+"\n"+str(para_unit_comb["unit procedure"].iloc[i])\
                        ,fontsize=10)
            ax.set_xlabel("")
            ax.set_ylabel("")

        fig.savefig(program+"QQplotScreening.pdf")
  


 #convert the map/filter/reduce into a list, for readability.
#list(map(....)) otherwise they are just map/filter objects which can not be displayed.
###### STR to NUM on a series
###### 1. str.isdigit(): a function only applied to string, and judge based on the string altogether. Does not judge if there is any digit contained in the string.
#eg: "100"-->T, "100+"-->F
###### 2. re.findall
#return all blocks of number

#re.findall('\d+', "00>=6000asdf90.sdf9")-->['00', '6000', '90', '9']

#return all sigle digits in a string(as a list)

#re.findall('\d', "00>=6000asdf90.sdf9")-->['0', '0', '6', '0', '0', '0', '9', '0', '9']

#Requires one number in each cell.
#A way to exam digit element-wise: use filter(treat a string not as a whole but as a list of elements)

#list(filter(str.isdigit, '200 grams'))  output:["2","0","0"]

#concat all elements in a list: ''.join(list(filter(str.isdigit, '200 grams')))

#convert this to a number: float(''.join(list(filter(str.isdigit, '200 grams'))))

#This function can be mapped to UAL/LAL for str->num conversion ONLY IF: 
    #There is only one number in each cell. There has to be at least one digit for "float" function to work.
#def string_to_num(input_string):
    
    #The input of this function has to be STRING; will raise error in filter otherwise
   # output=''.join(list(filter(str.isdigit, input_string)))
    #if (output==""): 
       # return(".")
    #else:
       # return(float(output))
    
#def string_to_num_re(input_string):
    #output_list=re.findall('\d+',input_string)
    #if (output_list==[]):
        #if the original cell/spec is missing/unspecified
       # return(".")
    #else:
        #take the only number in the list
        #return(float(output_list[0]))
        
        