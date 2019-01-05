import pandas as pd
import numpy as np

import os

from scipy.stats import anderson
from scipy.stats import kstest
from scipy.stats import shapiro

class APR_analytical_pipeline:
    def __init__(self....):
        
        
        
    def kurt(values):
        n=len(values)
        if (n<4): 
            result=np.nan
        else:
            sample_mean=np.mean(values)
            sample_std=np.std(values,ddof=1)
            result=np.sum(((values-sample_mean)/sample_std)**4)*n*(n+1)/((n-1)*(n-2)*(n-3))-3*(n-1)**2/((n-2)*(n-3))
        return(result)

    def skewness(values):
        n=len(values)
        if (n<3):
            result=np.nan
        else:
            sample_mean=np.mean(values)
            sample_std=np.std(values,ddof=1)
            result=n*np.sum(((values-sample_mean)/sample_std)**3)/((n-1)*(n-2))
        return(result)

    def CV(values):
        sample_mean=np.mean(values)
        sample_std=np.std(values,ddof=1)
        return(100*sample_std/sample_mean)

    def ppk_single_row(values,ual,lal,spec_status):
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

    def add_specside(df):
        #no spec default
        df["spec status"]=None
        #two specs
        df["spec status"].loc[(df.UAL!=".") & (df.LAL!=".")]=0
        #only lower spec
        df["spec status"].loc[(df.UAL==".") & (df.LAL!=".")]=-1
        #only upper spec
        df["spec status"].loc[(df.UAL!=".") & (df.LAL==".")]=1

    def add_clside(df):
        #no spec default
        df["ctrl status"]=None
        #two specs
        df["ctrl status"].loc[(df.LCL!=".") & (df.UCL!=".")]=0
        #only lower spec
        df["ctrl status"].loc[(df.UCL==".") & (df.LCL!=".")]=-1
        #only upper spec
        df["ctrl status"].loc[(df.UCL!=".") & (df.LCL==".")]=1

    def df_colname_unify(df):
        df=df.rename(columns={"ual":"UAL","lal":"LAL","lcl":"LCL","ucl":"UCL"})

    def quanti_quali_namelists(df):

        qualitative_paras=list(pd.unique(df["parameter name"][df['result type']=="TEXT"]))

        quantitative_paras=list(pd.unique(df["parameter name"][df['result type']=="Number"]))

        return(quantitative_paras,qualitative_paras)

    def para_unit_combi_df(df):
        para_unit_combi=df[["parameter name","unit procedure"]].copy().drop_duplicates().reset_index(drop=True)
        para_unit_combi=para_unit_combi.sort_values(by=["unit procedure","parameter name"])
        return(para_unit_combi)

    #only use a quantitative dataset as the input
    def spec_correction_rawset(df,para_unit_pair_df):
        for i in range(len(para_unit_pair_df)):
            para=para_unit_pair_df["parameter name"].iloc[i]
            unit=para_unit_pair_df["unit procedure"].iloc[i]

            sub_df_len=len(df.loc[(df['parameter name']==para) &(df["unit procedure"]==unit)])

            #more than one row, needs spec correction
            if (sub_df_len>1):

                #could loop using the uniq combination of unit and para, otherwise add this row number check to avoid void error
                ual=df[(df['parameter name']==para) &(df["unit procedure"]==unit)].UAL.iloc[-1]

                lal=df[(df['parameter name']==para) &(df["unit procedure"]==unit)].LAL.iloc[-1]
                df["UAL"].loc[(df['parameter name']==para) &(df["unit procedure"]==unit)]=ual
                df["LAL"].loc[(df['parameter name']==para) &(df["unit procedure"]==unit)]=lal

    def find_oos_index(df):
        df_twospecs=df[df["spec status"]==0]
        df_uspecs=df[df["spec status"]==1]
        df_lspecs=df[df["spec status"]==-1]

        #Find the indexes of OOS rows, try to do this in SQL
        twospecs_rows_oos=list(df_twospecs.index[(df_twospecs["value_num(recorded/full precision)"]<df_twospecs.LAL)\
                        |(df_twospecs["value_num(recorded/full precision)"]>df_twospecs.UAL)])

        uspecs_rows_oos=list(df_uspecs.index[(df_uspecs["value_num(recorded/full precision)"]>df_uspecs.UAL)])

        lspecs_rows_oos=list(df_lspecs.index[(df_lspecs["value_num(recorded/full precision)"]<df_lspecs.LAL)])

        return(list(set(twospecs_rows_oos+uspecs_rows_oos+lspecs_rows_oos)))

    def find_oot_index(df):
        df_twocls=df[df["ctrl status"]==0]
        df_ucl=df[df["ctrl status"]==1]
        df_lcl=df[df["ctrl status"]==-1]

        twocls_rows_oot=list(df_twocls.index[(df_twocls["value_num(recorded/full precision)"]<df_twocls.LCL)\
                    |(df_twocls["value_num(recorded/full precision)"]>df_twocls.UCL)])

        ucls_rows_oot=list(df_ucl.index[(df_ucl["value_num(recorded/full precision)"]>df_ucl.UCL)])

        lcls_rows_oot=list(df_lcl.index[(df_lcl["value_num(recorded/full precision)"]<df_lcl.LCL)])

        return(list(set(twocls_rows_oot+ucls_rows_oot+lcls_rows_oot)))

    #change the threshold(30) if we allow smaller sample size after excluding the OOS
    #input_type can be "" or "_no_OOS" or "no_OOST", etc, used to define the the column name
    def aggre_summary_df_generator(df, para_unit_pair_df, input_type):

        #initialize the summary output df
        aggre_stats_output_df=pd.DataFrame(columns=["parameter name",'unit procedure',#"PRDS",#\
                                               "N","mean","median","STD","variance","UAL","LAL","spec status","Ppk"+input_type,\
                                               "skewness","kurtosis","CV",#"anderson-darling"\
                                                "shapiro","KS"\
                                                  ])
        for i in range(len(para_unit_pair_df)):

            para=para_unit_pair_df["parameter name"].iloc[i]
            unit=para_unit_pair_df["unit procedure"].iloc[i]

            single_df=df.loc[(df['parameter name']==para) &(df["unit procedure"]==unit)].copy().sort_values(by="parameter date") 
            single_df_values=single_df["value_num(recorded/full precision)"]  

            ual=single_df.UAL.iloc[-1]
            lal=single_df.LAL.iloc[-1]
            spec_status=single_df["spec status"].iloc[-1]

            if (len(single_df)<30):
                ppk="NA1"

            if( (len(single_df)>=30) and (len(pd.unique(single_df_values)<5)) ):
                ppk="NA2"

            #NA3: no specs is addressed inside the ppk_single_row

            if ((len(single_df)>=30) and (len(pd.unique(single_df_values)>=5)) ):

                #prds=pd.unique(single_df["parameter detail"])

                ppk=ppk_single_row(single_df_values,ual,lal,spec_status)
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
                            "skewness":skewness(single_df_values),\
                            "kurtosis":kurt(single_df_values),\
                            "CV":CV(single_df_values),\
                            #"Anderson-Darling":anderson(sub_df),\
                            "shapiro":shap,\
                            #the p value of shapiro
                           "KS":ks},\
                            #the p value of kstest
                                        ignore_index=True)
        return(aggre_stats_output_df)

    def JMP_single_para_subset_generator(para_unit_pair_df, raw_df, folder_path):
        os.chdir(folder_path)
        for i in range(len(para_unit_pair_df)):
            para=para_unit_pair_df["parameter name"].iloc[i]
            unit=para_unit_pair_df["unit procedure"].iloc[i]
            single_subset=raw_df[(raw_df["parameter name"]==para)&(raw_df["unit procedure"]==unit)].copy()
            single_subset.to_csv(para+"-"+unit+"JMP.csv", index=False)

        return