#IVs#
#Body mass index ieu-a-2
#Height ieu-a-1032
#Waist-to-hip rati ieu-a-72
#Body fat ieu-a-999
#LDL cholesterol ieu-b-110
#HDL cholesterol ieu-b-109
#triglycerides ieu-b-111
#Total cholesterol ieu-a-301
#apolipoprotein A-I ieu-b-107
#Adiponectin ieu-a-1
#Fatty liver finn-b-NAFLD
#T2DM ebi-a-GCST006867
#HbA1c bbj-a-26
#Serum 25-Hydroxyvitamin D levels ebi-a-GCST90000618
#Uric acid bbj-a-57
#hypertension ukb-b-12493
#Systolic blood pressure ieu-b-38
#Diastolic blood pressure ieu-b-39

library(MRInstruments)
library(TwoSampleMR)
#outcome#
#finn-b-CD2_BENIGN_THYROID Benign neoplasm of thyroid gland#
###finn-b-C3_THYROID_GLAND  Malignant neoplasm of thyroid gland ncase989 ncontrol217,803#
exp_dat <-extract_instruments(outcomes=c('ieu-a-2',
                                         'ieu-a-1032',
                                         'ieu-a-72',
                                         'ieu-a-999',
                                         'ieu-b-110',
                                         'ieu-b-109',
                                         'ieu-b-111',
                                         'ieu-a-301',
                                         'ieu-b-107',
                                         'ieu-a-1',
                                         'finn-b-NAFLD',
                                         'ebi-a-GCST006867',
                                         'bbj-a-26',
                                         'ebi-a-GCST90000618',
                                         'bbj-a-57',
                                         'ukb-b-12493',
                                         'ieu-b-38',
                                         'ieu-b-39'))
outcome_data <- extract_outcome_data(snps = exp_dat$SNP, outcomes = "finn-b-C3_THYROID_GLAND",access_token = NULL)
dat <- harmonise_data(exp_dat, outcome_data, action = 2)


res <- mr(dat)
res2 <- subset_on_method(res) 
#default isto subset on either the IVW method (>1 instrumental SNP) or Wald ratiomethod (1 instrumental SNP).
res3 <- sort_1_to_many(res,b="b",sort_action=4)
#this sorts results by decreasing effect size (largest effect at top of theplot)
res4 <- split_exposure(res) 
# to keep the Yaxis label clean we exclude the exposure ID labels from the exposure column
res$weight <- 1/res$se
res$low<-res$b-1.96*res$se
res$up<-res$b+1.96*res$se
res$low95<-exp(res$b-1.96*res$se)
res$high95<-exp(res$b+1.96*res$se)

min(exp(res$b-1.96*res$se)) 

#identify value for 'lo' in forest_plot_1_to_many
max(exp(res$b+1.96*res$se))
generate_odds_ratios(res)
#identify valuefor 'up' in forest_plot_1_to_many
forest_plot_1_to_many(res,b="b",se="se",
                      exponentiate=T,ao_slc=F,lo=0.3,up=2.5,
                      TraitM="exposure",col1_width=2,by=NULL, 
                      trans="log2",xlab="OR for CHD per SD increase in risk factor (95% confidenceinterval)",weight="weight")

res$pval <- formatC(res$pval, format ="e", digits = 2)

p1<-forest_plot_1_to_many(res,b="b",se="se",
                          exponentiate=T,ao_slc=F,lo=0.3,up=2.5,
                          TraitM="exposure",by=NULL,
                          trans="log2",xlab="OR for CHD per SD increase in riskfactor (95% CI)",   
                          weight="weight",subheading_size=11,
                          col1_title="Risk factor",
                          col1_width=2.5,
                          col_text_size=4,
                          addcols=c("nsnp","pval"),
                          addcol_widths=c(1.0,1.0)
)
res<-mr(dat)
res<-split_exposure(res) 
# to keep the Yaxis label clean we exclude the exposure ID labels from the exposure column
res<-sort_1_to_many(res,group="exposure",sort_action=3,priority="Inverse variance weighted",trait_m="method")
p<-forest_plot_1_to_many(res,b="b",se="se",
                         exponentiate=T,trans="log2",ao_slc=F,lo=0.03,
                         up=22,col1_width=4,by="exposure",TraitM="method",
                         # xlab=######
                         subheading_size=12,col_text_size=4,
                         addcols=c("nsnp","pval"),
                         addcol_widths=c(1.0,1.0))

library(ggplot2)

exp_dat <-extract_instruments(outcomes="ieu-b-39")

outcome_data <- extract_outcome_data(snps = exp_dat$SNP, outcomes = "finn-b-C3_THYROID_GLAND",access_token = NULL)

H_data <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = outcome_data 
)

mr_results<-mr(H_data)
mr_results
mr(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
mr_method_list()
head(mr_method_list())[,1:2]
generate_odds_ratios(mr_results)
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
plot1 <- mr_scatter_plot(mr_results, H_data)
plot1
res_single <- mr_singlesnp(H_data)
plot2 <- mr_forest_plot(res_single)
plot2
res_loo <- mr_leaveoneout(H_data)
plot3 <- mr_leaveoneout_plot(res_loo)
plot3
plot4 <- mr_funnel_plot(res_single)
plot4


#multi-MR
##ieu-b-110
#HDL cholesterol ieu-b-109
#triglycerides ieu-b-111
#Total cholesterol ieu-a-301###
id_exposure <- c('ieu-b-110','ieu-b-109','ieu-b-111','ieu-b-39')
id_outcome <- "finn-b-C3_THYROID_GLAND"
exposure_dat <- mv_extract_exposures(id_exposure)
dim(exposure_dat)
oucome_dat <- extract_outcome_data(exposure_dat$SNP, id_outcome) 

mvdat <- mv_harmonise_data(exposure_dat,  oucome_dat) 
res_mult <- mv_multiple(mvdat) 
res_mult 

##OR and 95%CI
res_mult$result$weight <- 1/res_mult$result$se
res_mult$result$low<-res_mult$result$b-1.96*res_mult$result$se
res_mult$result$up<-res_mult$result$b+1.96*res_mult$result$se
res_mult$result$low95<-exp(res_mult$result$b-1.96*res_mult$result$se)
res_mult$result$high95<-exp(res_mult$result$b+1.96*res_mult$result$se)
res_mult$result$or<-exp(res_mult$result$b)


