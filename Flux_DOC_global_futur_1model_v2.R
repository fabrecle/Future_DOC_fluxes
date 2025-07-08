# Load packages ####
library(plyr)
library(broom)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(matrixStats)
library(Kendall)
library(trend)
library(sf)
library(rworldmap)
# library(raster)
library(spData)
library(tmap)
library(nlmrt)
library(nls2)
library(envalysis)
library(patchwork)
library(foreach)
library(doParallel)
library(paletteer)
library("gridExtra")

# Directory ####
dir = 'E:/Papiers/2020_Fabre et al Global carbon futur/GRDC'
# Preparing dataset ####
climats = as.data.frame(read_xlsx("E:/Papiers/2020_Fabre et al Global carbon futur/Climats_future_v3.xlsx", sheet="Climats_future"))

regions = as.data.frame(read_xlsx("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/regions.xlsx"))
regions = regions[regions$SUM==45,]

climats = merge(climats, regions[c('region','wribasin_u')], by.x = 'NAME', by.y="wribasin_u")

rownames(climats) <- climats[,1]
climats[,1] <- NULL
climats$OBJECTID <- NULL

climats = na.omit(climats)

variables = c('Q', 'PCP', 'TMP', 'SOC')

periods = c('hist', 'fut')

models = c('gfdl', 'hadgem', 'ipsl', 'miroc', 'noresm')
RCPs = c('rcp26', 'rcp85')

type = c('mean', 'sd')


all_climats = c('A','B','C','D','E')
all_climats_long = c('Tropical', 'Semi-Arid', 'Temperate', 'Cold', 'Polar')
all_oceans = c('Atlantic Ocean','Arctic Ocean', 'Indian Ocean', 'Pacific Ocean')
all_oceans_short = c('Atl','Arc','Ind','Pac')

new_outputs = expand.grid('DOC',periods,RCPs,models,type)
new_outputs = new_outputs[order(new_outputs$Var1,new_outputs$Var2, new_outputs$Var3, new_outputs$Var4),]
new_variables = expand.grid(variables,periods,RCPs,models)
new_variables = new_variables[order(new_variables$Var1,new_variables$Var2, new_variables$Var3),]

climats[,paste(new_outputs$Var1,new_outputs$Var2,new_outputs$Var3,new_outputs$Var4, new_outputs$Var5, sep='_')] = NA
climats[,paste(new_variables$Var1,new_variables$Var2,new_variables$Var3,new_variables$Var4, sep='_')] = NA


path = "E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Files"
path_final = "E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final"

path_files = "E:/Papiers/2020_Fabre et al Global carbon futur/Future_discharge/Complete_Routing_run_"


models_Q = c('gfdl-esm2m', 'hadgem2-es', 'ipsl-cm5a-lr', 'miroc-esm-chem', 'noresm1-m')
RCPs_Q = c('HistAndRcp2p6', 'HistAndRcp8p5')


# Regression coefficients correction Monte Carlo ####
data_alpha = read_xlsx("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/alpha_SOC_PCP_GLM.xlsx")

MC_reg = function(data, y, sd_y, SOC, PRECIP, N_MC){
  a=c(); b=c(); r=c(); pval=c()
  for(i in 1:N_MC){
    alpha_MC = rnorm(nrow(data), mean=data[[y]], sd=data[[sd_y]])
    entry_MC = rnorm(nrow(data), mean=data[[SOC]]/(data[[PRECIP]]/12),
                     sd=(data[[SOC]]/(data[[PRECIP]]/12))*(0.1 + 0.1))

    P = lm(alpha_MC~entry_MC)
    a[i] = summary(P)$coefficients[2]
    b[i] = summary(P)$coefficients[1]
    r[i] = glance(P)$r.squared
    pval[i] = glance(P)$p.value
  }
  return(data.frame(a=a, b=b, r=r, pval=pval))
}

GFDL_MC = MC_reg(data_alpha, 'alpha', 'std_alpha', 'GFDL_SOC', 'GFDL_PCP', 10000)
Hadgem_MC = MC_reg(data_alpha, 'alpha', 'std_alpha', 'Hadgem_SOC', 'Hadgem_PCP', 10000)
IPSL_MC = MC_reg(data_alpha, 'alpha', 'std_alpha', 'IPSL_SOC', 'IPSL_PCP', 10000)
MIROC_MC = MC_reg(data_alpha, 'alpha', 'std_alpha', 'MIROC_SOC', 'MIROC_PCP', 10000)
NorESM_MC = MC_reg(data_alpha, 'alpha', 'std_alpha', 'NorESM_SOC', 'NorESM_PCP', 10000)

summary_a_b = data.frame(slope=c(mean(GFDL_MC$a),mean(Hadgem_MC$a),mean(IPSL_MC$a),mean(MIROC_MC$a),mean(NorESM_MC$a)),
                    slope_std=c(sd(GFDL_MC$a),sd(Hadgem_MC$a),sd(IPSL_MC$a),sd(MIROC_MC$a),sd(NorESM_MC$a)),
                    intercept=c(mean(GFDL_MC$b),mean(Hadgem_MC$b),mean(IPSL_MC$b),mean(MIROC_MC$b),mean(NorESM_MC$b)),
                    intercept_std=c(sd(GFDL_MC$b),sd(Hadgem_MC$b),sd(IPSL_MC$b),sd(MIROC_MC$b),sd(NorESM_MC$b)),
                    R = c(mean(GFDL_MC$r),mean(Hadgem_MC$r),mean(IPSL_MC$r),mean(MIROC_MC$r),mean(NorESM_MC$r)),
                    R_std = c(sd(GFDL_MC$r),sd(Hadgem_MC$r),sd(IPSL_MC$r),sd(MIROC_MC$r),sd(NorESM_MC$r)),
                    pval=c(mean(GFDL_MC$pval),mean(Hadgem_MC$pval),mean(IPSL_MC$pval),mean(MIROC_MC$pval),mean(NorESM_MC$pval)),
                    row.names = c('gfdl','hadgem','ipsl','miroc','noresm'))



MC_plot = function(data, model, model_title, MC){
  p = ggplot(data, aes(y=.data[['alpha']], x=.data[[paste0(model,'_SOC')]]/(.data[[paste0(model,'_PCP')]]/12))) +
    xlab(expression(frac("SOC content (kg."~m^-3*")", "Precipitation (mm."~m^-1*")"))) +
    ylab(expression(alpha)) +
    ggtitle(paste(model_title, " (R²=",round(summary_a_b[tolower(model),]$R,2), ifelse(summary_a_b[tolower(model),]$pval<0.05,", p<0.05",paste0(', p=',round(summary_a_b[tolower(model),]$pval,2))),")")) +
    geom_point(aes(y=.data[['alpha']]), na.rm = FALSE) +
    expand_limits(x=0, y=0) +
    geom_line(data=data, aes(x=.data[[paste0(model,'_SOC')]]/(.data[[paste0(model,'_PCP')]]/12), y=mean(MC$a)*.data[[paste0(model,'_SOC')]]/(.data[[paste0(model,'_PCP')]]/12)+mean(MC$b)), linewidth=1) +
    geom_ribbon(data=data, aes(x=.data[[paste0(model,'_SOC')]]/(.data[[paste0(model,'_PCP')]]/12), ymin=(mean(MC$a)-sd(MC$a))*.data[[paste0(model,'_SOC')]]/(.data[[paste0(model,'_PCP')]]/12)+(mean(MC$b)-sd(MC$b)),
                                     ymax=(mean(MC$a)+sd(MC$a))*.data[[paste0(model,'_SOC')]]/(.data[[paste0(model,'_PCP')]]/12)+(mean(MC$b)+sd(MC$b))), inherit.aes = FALSE, fill="darkgrey", alpha=0.3) +
    theme_publish(base_size = 10, base_linewidth = 0.7)

  p_a = ggplot(MC, aes(x=a)) +
    ylab("Frequency") +
    geom_freqpoly() +
    geom_vline(xintercept=mean(MC$a) + sd(MC$a), linetype="dashed") +
    geom_vline(xintercept=mean(MC$a) - sd(MC$a), linetype="dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_publish(base_size = 10, base_linewidth = 0.7)

  p_b = ggplot(MC, aes(x=b)) +
    ylab("Frequency") +
    geom_freqpoly() +
    geom_vline(xintercept=mean(MC$b) + sd(MC$b), linetype="dashed") +
    geom_vline(xintercept=mean(MC$b) - sd(MC$b), linetype="dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_publish(base_size = 10, base_linewidth = 0.7)

  return(p / (p_a | p_b) + plot_layout(nrow = 2, heights = c(3, 1)))
}

ggarrange(MC_plot(data_alpha, 'GFDL', 'GFDL-ESM2M', GFDL_MC),
          MC_plot(data_alpha, 'Hadgem', 'HadGEM2-ES', Hadgem_MC),
          MC_plot(data_alpha, 'IPSL', 'IPSL-CM5A-LR', IPSL_MC),
          MC_plot(data_alpha, 'MIROC', 'MIROC-ESM-CHEM', MIROC_MC),
          MC_plot(data_alpha, 'NorESM', 'NorESM1-M', NorESM_MC), ncol = 2, nrow=3)



# Creating DOC files per watershed : MEAN, MIN and MAX ####
# sum = 0; area_tot = 0; iter = 0
# 
# columns_year = seq(1971,2099,1)
# columns_julian = seq(1,366,1)
# 
# 
# # Yearly initialization
# climats_models = expand.grid(RCPs,all_climats,models,type)
# climats_models = climats_models[order(climats_models$Var1, climats_models$Var2, climats_models$Var3),]
# oceans_models = expand.grid(RCPs,all_oceans_short,models,type)
# oceans_models = oceans_models[order(oceans_models$Var1, oceans_models$Var2, oceans_models$Var3),]
# 
# results_climats = data.frame(year=seq(1971,2099,1))
# results_climats[,paste(climats_models$Var1,climats_models$Var2,climats_models$Var3,climats_models$Var4,sep='_')] = double(129)
# 
# results_oceans = data.frame(year=seq(1971,2099,1))
# results_oceans[,paste(oceans_models$Var1,oceans_models$Var2,oceans_models$Var3,oceans_models$Var4,sep='_')] = double(129)
# 
# df_models = data.frame(Y = seq(1971,2099,1))
# 
# # Yearly discharge
# results_climats_Q = data.frame(year=seq(1971,2099,1))
# results_climats_Q[,paste(climats_models$Var1,climats_models$Var2,climats_models$Var3,climats_models$Var4,sep='_')] = double(129)
# 
# results_oceans_Q = data.frame(year=seq(1971,2099,1))
# results_oceans_Q[,paste(oceans_models$Var1,oceans_models$Var2,oceans_models$Var3,oceans_models$Var4,sep='_')] = double(129)
# 
# df_models_Q = data.frame(Y = seq(1971,2099,1))
# 
# 
# # Daily initialization
# daily_climats_hist = data.frame(Julian=seq(1,366, by=1))
# daily_climats_hist[,paste(climats_models$Var1,climats_models$Var2,climats_models$Var3,climats_models$Var4,sep='_')] = double(366)
# daily_oceans_hist = data.frame(Julian=seq(1,366, by=1))
# daily_oceans_hist[,paste(oceans_models$Var1,oceans_models$Var2,oceans_models$Var3,oceans_models$Var4,sep='_')] = double(366)
# daily_models_hist = data.frame(Julian=seq(1,366, by=1))
# 
# daily_climats = data.frame(Julian=seq(1,366, by=1))
# daily_climats[,paste(climats_models$Var1,climats_models$Var2,climats_models$Var3,climats_models$Var4,sep='_')] = double(366)
# daily_oceans = data.frame(Julian=seq(1,366, by=1))
# daily_oceans[,paste(oceans_models$Var1,oceans_models$Var2,oceans_models$Var3,oceans_models$Var4,sep='_')] = double(366)
# daily_models = data.frame(Julian=seq(1,366, by=1))
# 
# 
# MC_DOCflux = function(Q, area, summary, model, N_MC){
# 
#   cl = makeCluster(parallel::detectCores() - 1)
#   registerDoParallel(cl)
# 
#   alpha_MC = rnorm(N_MC, mean = summary[model,3], sd = summary[model,4]) +
#               rnorm(N_MC, mean = summary[model,1], sd = summary[model,2])*
#               (apply(soilc_sel, 2, FUN = function(x) rnorm(N_MC, mean = x, sd = x/10))/
#                  (apply(pcp_sel, 2, FUN = function(x) rnorm(N_MC, mean = x, sd = x/10))/12))
# 
#   tmp_sel_MC = apply(tmp_sel, 2, FUN = function(x) rnorm(N_MC, mean = x, sd = 1))
#   beta_MC = ifelse(tmp_sel_MC<(-13), 4, (12.6 / (tmp_sel_MC + 16.1)) - 0.03)
# 
#   oper = foreach(l=1:N_MC, .packages=c('plyr','dplyr')) %dopar% {
#     Q_MC = data.frame(Y=Q$Y, julian=Q$julian, Value = rnorm(nrow(Q), mean=Q$Value, sd=Q$Value/10))
#     alpha_df = data.frame(Y = 1971:2099, alpha = alpha_MC[l,])
#     beta_df = data.frame(Y = 1971:2099, beta = beta_MC[l,])
# 
#     df_calc = Reduce(function(x, y) merge(x, y, by='Y', all=TRUE), list(Q_MC, alpha_df, beta_df))
#     df_calc$Value_j = df_calc$Value*86400*1000
#     df_calc$DOC_flux = df_calc$alpha * (df_calc$Value_j/area) / (df_calc$beta + (df_calc$Value_j/area)) *1e-15* df_calc$Value_j
# 
#     result1 = ddply(df_calc,~Y,summarise, FDOC = sum(DOC_flux, na.rm=T))$FDOC
#     result2 = ddply(df_calc[between(df_calc$Y,1971,2005),],~julian,summarise, FDOC = mean(DOC_flux, na.rm=T))$FDOC
#     result3 = ddply(df_calc[between(df_calc$Y,2071,2099),],~julian,summarise, FDOC = mean(DOC_flux, na.rm=T))$FDOC
# 
#     return(list(result1,result2,result3))
#   }
# 
#   stopCluster(cl)
# 
#   result1 = as.data.frame(do.call(rbind,lapply(oper,function(x){x[[1]]})))
#   result2 = as.data.frame(do.call(rbind,lapply(oper,function(x){x[[2]]})))
#   result3 = as.data.frame(do.call(rbind,lapply(oper,function(x){x[[3]]})))
# 
#   result1_mean = colMeans(as.matrix(result1))
#   result2_mean = colMeans(as.matrix(result2))
#   result3_mean = colMeans(as.matrix(result3))
# 
#   result1_sd = colSds(as.matrix(result1))
#   result2_sd = colSds(as.matrix(result2))
#   result3_sd = colSds(as.matrix(result3))
# 
#   names(result1_mean) = columns_year
#   names(result2_mean) = columns_julian
#   names(result3_mean) = columns_julian
#   names(result1_sd) = columns_year
#   names(result2_sd) = columns_julian
#   names(result3_sd) = columns_julian
# 
#   return(list(result1_mean, result1_sd,
#               result2_mean, result2_sd,
#               result3_mean, result3_sd))
# }
# 
# 
# for (name in row.names(climats)){
#   data_watershed = climats[row.names(climats)==name,]
#   print(c(name,data_watershed$CLIMAT,data_watershed$OCEAN))
#   area = data_watershed$Area * 1e6
# 
#   # HISTORICAL
#   for (j in 1:length(RCPs)){
#     for (k in 1:length(models)){
# 
#       print(c(models[k],RCPs[j]))
# 
#       Q = read.table(file.path(paste0(path_files,models_Q[k],"_",RCPs_Q[j],"_1961-2099_v1_25bands"),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
#       Q = Q[Q$Y >= 1971,]
#       Q$julian = format(as.Date(with(Q, paste(Y, M, D,sep="-"))), "%j")
# 
#       Q_yearly = ddply(Q,~Y,summarise, Q_yearly = sum(Value*86400/1e9, na.rm=T))
# 
#       pcp_hist = read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/Hist/PCP_",models[k],".csv"))
#       pcp_fut = read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/PCP_",models[k],"_",RCPs[j],".csv"))
#       pcp = merge(pcp_hist, pcp_fut, by="region")
# 
#       pcp_sel = as.data.frame(pcp[pcp$region==data_watershed$region,])
#       pcp_sel$region = NULL; pcp_sel$X.x = NULL; pcp_sel$X.y = NULL; pcp_sel$X2100 = NULL
#       colnames(pcp_sel) = columns_year
# 
#       tmp_hist = read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/Hist/TMP_",models[k],".csv"))
#       tmp_fut = read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/TMP_",models[k],"_",RCPs[j],".csv"))
#       tmp = merge(tmp_hist, tmp_fut, by="region")
# 
#       tmp_sel = as.data.frame(tmp[tmp$region==data_watershed$region,])
#       tmp_sel$region = NULL; tmp_sel$X.x = NULL; tmp_sel$X.y = NULL; tmp_sel$X2100 = NULL
#       colnames(tmp_sel) = columns_year
# 
#       soilc_hist = read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/Hist/SOC_",models[k],".csv"))
#       soilc_fut = read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/SOC_",models[k],"_",RCPs[j],".csv"))
#       soilc = merge(soilc_hist, soilc_fut, by="region")
# 
#       if (models[k] %in% c('ipsl', 'miroc', 'noresm')){
#         litter_hist =  read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/Hist/Litter_",models[k],".csv"))
#         litter_fut = read.csv(paste0("E:/Papiers/2020_Fabre et al Global carbon futur/DL_variables_futur/Litter_",models[k],"_",RCPs[j],".csv"))
#         litter = merge(litter_hist, litter_fut, by="region")
#         soilc[3:length(soilc)] = soilc[3:length(soilc)] + litter[3:length(litter)]
#       }
# 
#       soilc_sel = as.data.frame(soilc[soilc$region==data_watershed$region,]) #kg/m²
#       soilc_sel$region = NULL; soilc_sel$X.x = NULL; soilc_sel$X.y = NULL; soilc_sel$X2100 = NULL
#       colnames(soilc_sel) = columns_year
# 
#       # Pour chaque variable et chaque année, tirer 1000 valeurs entre la valeur moyenne et ses STD
# 
#       DOCflux_MC = MC_DOCflux(Q, area, summary_a_b, k, 988)
#       
#       # save lists as csv ?
#       write.csv(data.frame(mean=DOCflux_MC[[1]],sd=DOCflux_MC[[2]]), file = paste0(path_final,"/",name,"_",models[k],"_",RCPs[j],".csv"))
# 
#       # fill climats for summary
#       climats[row.names(climats)==name,][[paste0("DOC_","hist_",RCPs[j],"_",models[k],'_mean')]] = mean(DOCflux_MC[[1]][names(DOCflux_MC[[1]])<=2000])
#       climats[row.names(climats)==name,][[paste0("DOC_","hist_",RCPs[j],"_",models[k],'_sd')]] = mean(DOCflux_MC[[2]][names(DOCflux_MC[[2]])<=2000])
#       climats[row.names(climats)==name,][[paste0("DOC_","fut_",RCPs[j],"_",models[k],'_mean')]] = mean(DOCflux_MC[[1]][names(DOCflux_MC[[1]])>=2071])
#       climats[row.names(climats)==name,][[paste0("DOC_","fut_",RCPs[j],"_",models[k],'_sd')]] = mean(DOCflux_MC[[2]][names(DOCflux_MC[[2]])>=2071])
# 
#       climats[row.names(climats)==name,][[paste0("Q_","hist_",RCPs[j],"_",models[k])]] = mean(Q_yearly[Q_yearly$Y<=2000,]$Q_year)
#       climats[row.names(climats)==name,][[paste0("TMP_","hist_",RCPs[j],"_",models[k])]] = mean(data.frame(t(tmp_sel))[as.integer(colnames(tmp_sel))<=2000,])
#       climats[row.names(climats)==name,][[paste0("PCP_","hist_",RCPs[j],"_",models[k])]] = mean(data.frame(t(pcp_sel))[as.integer(colnames(pcp_sel))<=2000,])
#       climats[row.names(climats)==name,][[paste0("SOC_","hist_",RCPs[j],"_",models[k])]] = mean(data.frame(t(soilc_sel))[as.integer(colnames(soilc_sel))<=2000,])
# 
#       climats[row.names(climats)==name,][[paste0("Q_","fut_",RCPs[j],"_",models[k])]] = mean(Q_yearly[Q_yearly$Y>=2071,]$Q_year)
#       climats[row.names(climats)==name,][[paste0("TMP_","fut_",RCPs[j],"_",models[k])]] = mean(data.frame(t(tmp_sel))[as.integer(colnames(tmp_sel))>=2071,])
#       climats[row.names(climats)==name,][[paste0("PCP_","fut_",RCPs[j],"_",models[k])]] = mean(data.frame(t(pcp_sel))[as.integer(colnames(pcp_sel))>=2071,])
#       climats[row.names(climats)==name,][[paste0("SOC_","fut_",RCPs[j],"_",models[k])]] = mean(data.frame(t(soilc_sel))[as.integer(colnames(soilc_sel))>=2071,])
# 
#       # fill climate and ocean columns
#       results_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] = results_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] + as.vector(DOCflux_MC[[1]])
#       results_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] = results_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] + as.vector(DOCflux_MC[[2]])
#       results_climats_Q[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] = results_climats_Q[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] + Q_yearly$Q_year
#       results_climats_Q[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] = results_climats_Q[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] + (Q_yearly$Q_year/10)
# 
#       daily_climats_hist[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] = daily_climats_hist[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] + as.vector(DOCflux_MC[[3]])
#       daily_climats_hist[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] = daily_climats_hist[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] + as.vector(DOCflux_MC[[4]])
#       daily_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] = daily_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_mean')]] + as.vector(DOCflux_MC[[5]])
#       daily_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] = daily_climats[[paste0(RCPs[j],"_",data_watershed$CLIMAT,"_",models[k],'_sd')]] + as.vector(DOCflux_MC[[6]])
# 
#       results_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] = results_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] + as.vector(DOCflux_MC[[1]])
#       results_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] = results_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] + as.vector(DOCflux_MC[[2]])
#       results_oceans_Q[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] = results_oceans_Q[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] + Q_yearly$Q_year
#       results_oceans_Q[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] = results_oceans_Q[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] + (Q_yearly$Q_year/10)
# 
#       daily_oceans_hist[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] = daily_oceans_hist[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] + as.vector(DOCflux_MC[[3]])
#       daily_oceans_hist[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] = daily_oceans_hist[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] + as.vector(DOCflux_MC[[4]])
#       daily_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] = daily_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_mean')]] + as.vector(DOCflux_MC[[5]])
#       daily_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] = daily_oceans[[paste0(RCPs[j],"_",substr(data_watershed$OCEAN, 1, 3),"_",models[k],'_sd')]] + as.vector(DOCflux_MC[[6]])
# 
#       }
#   }
# }
# 
# write.csv(climats, file=file.path(path_final,paste0("climats_new.csv")))
# write.csv(results_climats, file=file.path(path_final,paste0("results_climats_new.csv")))
# write.csv(results_climats_Q, file=file.path(path_final,paste0("results_climats_Q_new.csv")))
# write.csv(daily_climats_hist, file=file.path(path_final,paste0("daily_climats_hist_new.csv")))
# write.csv(daily_climats, file=file.path(path_final,paste0("daily_climats_new.csv")))
# write.csv(results_oceans, file=file.path(path_final,paste0("results_oceans_new.csv")))
# write.csv(results_oceans_Q, file=file.path(path_final,paste0("results_oceans_Q_new.csv")))
# write.csv(daily_oceans_hist, file=file.path(path_final,paste0("daily_oceans_hist_new.csv")))
# write.csv(daily_oceans, file=file.path(path_final,paste0("daily_oceans_new.csv")))


# Data_analysis ####

# Reimport datasets
climats = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/climats_new.csv')
rownames(climats) <- climats[,1]
climats[,1] <- NULL
results_climats = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/results_climats_new.csv')
results_climats_Q = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/results_climats_Q_new.csv')
daily_climats_hist = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/daily_climats_hist_new.csv')
daily_climats = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/daily_climats_new.csv')
results_oceans = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/results_oceans_new.csv')
results_oceans_Q = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/results_oceans_Q_new.csv')
daily_oceans_hist = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/daily_oceans_hist_new.csv')
daily_oceans = read.csv('E:/Papiers/2020_Fabre et al Global carbon futur/GRDC/Final/daily_oceans_new.csv')

# Results ####
results_climats_Q_av_sd = data.frame(year = results_climats_Q$year,
                                   rcp26_A_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_A.*mean$")))),
                                   rcp26_A_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_A.*sd$")))^2)),
                                   rcp26_B_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_B.*mean$")))),
                                   rcp26_B_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_B.*sd$")))^2)),
                                   rcp26_C_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_C.*mean$")))),
                                   rcp26_C_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_C.*sd$")))^2)),
                                   rcp26_D_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_D.*mean$")))),
                                   rcp26_D_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_D.*sd$")))^2)),
                                   rcp26_E_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_E.*mean$")))),
                                   rcp26_E_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp26_E.*sd$")))^2)),
                                   rcp85_A_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_A.*mean$")))),
                                   rcp85_A_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_A.*sd$")))^2)),
                                   rcp85_B_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_B.*mean$")))),
                                   rcp85_B_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_B.*sd$")))^2)),
                                   rcp85_C_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_C.*mean$")))),
                                   rcp85_C_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_C.*sd$")))^2)),
                                   rcp85_D_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_D.*mean$")))),
                                   rcp85_D_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_D.*sd$")))^2)),
                                   rcp85_E_mean = rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_E.*mean$")))),
                                   rcp85_E_sd = sqrt(rowMeans(as.data.frame(results_climats_Q %>% select(matches("^rcp85_E.*sd$")))^2)))


results_all_Q = data.frame(year = results_climats_Q_av_sd$year,
                         rcp26 = rowSums(as.data.frame(results_climats_Q_av_sd %>% select(matches("^rcp26.*mean$")))),
                         rcp26_sd = rowSums(as.data.frame(results_climats_Q_av_sd %>% select(matches("^rcp26.*sd$")))),
                         rcp85 = rowSums(as.data.frame(results_climats_Q_av_sd %>% select(matches("^rcp85.*mean$")))),
                         rcp85_sd = rowSums(as.data.frame(results_climats_Q_av_sd %>% select(matches("^rcp85.*sd$")))))


data.frame(control_global = mean(mean(results_all_Q[results_all_Q$year <=2000,]$rcp26),mean(results_all_Q[results_all_Q$year <=2000,]$rcp85)),
           control_global_sd = sqrt(mean(mean(results_all_Q[results_all_Q$year <= 2000,]$rcp26_sd^2),
                                         mean(results_all_Q[results_all_Q$year <= 2000,]$rcp85_sd^2))),
           future_rcp26 = mean(results_all_Q[results_all_Q$year >=2071,]$rcp26),
           future_rcp26_sd = sqrt(mean(results_all_Q[results_all_Q$year >=2071,]$rcp26_sd^2)),
           future_rcp85 = mean(results_all_Q[results_all_Q$year >=2071,]$rcp85),
           future_rcp85_sd = sqrt(mean(results_all_Q[results_all_Q$year >=2071,]$rcp85_sd^2)))

data.frame(control_A = mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp26_A_mean),mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp85_A_mean)),
           control_A_sd = sqrt(mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp26_A_sd^2),
                                         mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp85_A_sd^2))),
           control_B = mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp26_B_mean),mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp85_B_mean)),
           control_B_sd = sqrt(mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp26_B_sd^2),
                                    mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp85_B_sd^2))),
           control_C = mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp26_C_mean),mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp85_C_mean)),
           control_C_sd = sqrt(mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp26_C_sd^2),
                                    mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp85_C_sd^2))),
           control_D = mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp26_D_mean),mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp85_D_mean)),
           control_D_sd = sqrt(mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp26_D_sd^2),
                                    mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp85_D_sd^2))),
           control_E = mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp26_E_mean),mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <=2000,]$rcp85_E_mean)),
           control_E_sd = sqrt(mean(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp26_E_sd^2),
                                    mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year <= 2000,]$rcp85_E_sd^2))),
           
           future_A_rcp26 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_A_mean),
           future_A_rcp26_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_A_sd^2)),
           future_B_rcp26 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_B_mean),
           future_B_rcp26_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_B_sd^2)),
           future_C_rcp26 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_C_mean),
           future_C_rcp26_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_C_sd^2)),
           future_D_rcp26 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_D_mean),
           future_D_rcp26_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_D_sd^2)),
           future_E_rcp26 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_E_mean),
           future_E_rcp26_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp26_E_sd^2)),
           
           future_A_rcp85 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_A_mean),
           future_A_rcp85_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_A_sd^2)),
           future_B_rcp85 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_B_mean),
           future_B_rcp85_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_B_sd^2)),
           future_C_rcp85 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_C_mean),
           future_C_rcp85_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_C_sd^2)),
           future_D_rcp85 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_D_mean),
           future_D_rcp85_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_D_sd^2)),
           future_E_rcp85 = mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_E_mean),
           future_E_rcp85_sd = sqrt(mean(results_climats_Q_av_sd[results_climats_Q_av_sd$year >=2071,]$rcp85_E_sd^2))
           )

results_oceans_Q_av_sd = data.frame(year = results_oceans_Q$year,
                                     rcp26_Atl_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Atl.*mean$")))),
                                     rcp26_Atl_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Atl.*sd$")))^2)),
                                     rcp26_Arc_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Arc.*mean$")))),
                                     rcp26_Arc_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Arc.*sd$")))^2)),
                                     rcp26_Ind_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Ind.*mean$")))),
                                     rcp26_Ind_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Ind.*sd$")))^2)),
                                     rcp26_Pac_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Pac.*mean$")))),
                                     rcp26_Pac_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp26_Pac.*sd$")))^2)),
                                     rcp85_Atl_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Atl.*mean$")))),
                                     rcp85_Atl_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Atl.*sd$")))^2)),
                                     rcp85_Arc_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Arc.*mean$")))),
                                     rcp85_Arc_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Arc.*sd$")))^2)),
                                     rcp85_Ind_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Ind.*mean$")))),
                                     rcp85_Ind_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Ind.*sd$")))^2)),
                                     rcp85_Pac_mean = rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Pac.*mean$")))),
                                     rcp85_Pac_sd = sqrt(rowMeans(as.data.frame(results_oceans_Q %>% select(matches("^rcp85_Pac.*sd$")))^2)))

data.frame(control_Atl = mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp26_Atl_mean),mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp85_Atl_mean)),
           control_Atl_sd = sqrt(mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp26_Atl_sd^2),
                                    mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp85_Atl_sd^2))),
           control_Arc = mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp26_Arc_mean),mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp85_Arc_mean)),
           control_Arc_sd = sqrt(mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp26_Arc_sd^2),
                                    mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp85_Arc_sd^2))),
           control_Ind = mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp26_Ind_mean),mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp85_Ind_mean)),
           control_Ind_sd = sqrt(mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp26_Ind_sd^2),
                                    mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp85_Ind_sd^2))),
           control_Pac = mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp26_Pac_mean),mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <=2000,]$rcp85_Pac_mean)),
           control_Pac_sd = sqrt(mean(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp26_Pac_sd^2),
                                    mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year <= 2000,]$rcp85_Pac_sd^2))),
           
           future_Atl_rcp26 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Atl_mean),
           future_Atl_rcp26_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Atl_sd^2)),
           future_Arc_rcp26 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Arc_mean),
           future_Arc_rcp26_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Arc_sd^2)),
           future_Ind_rcp26 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Ind_mean),
           future_Ind_rcp26_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Ind_sd^2)),
           future_Pac_rcp26 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Pac_mean),
           future_Pac_rcp26_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp26_Pac_sd^2)),
           
           future_Atl_rcp85 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Atl_mean),
           future_Atl_rcp85_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Atl_sd^2)),
           future_Arc_rcp85 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Arc_mean),
           future_Arc_rcp85_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Arc_sd^2)),
           future_Ind_rcp85 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Ind_mean),
           future_Ind_rcp85_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Ind_sd^2)),
           future_Pac_rcp85 = mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Pac_mean),
           future_Pac_rcp85_sd = sqrt(mean(results_oceans_Q_av_sd[results_oceans_Q_av_sd$year >=2071,]$rcp85_Pac_sd^2))
)

# dataframe by ocean
results_oceans_av_sd = data.frame(year = results_oceans$year,
                                    rcp26_Atl_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Atl.*mean$")))),
                                    rcp26_Atl_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Atl.*sd$")))^2)),
                                    rcp26_Arc_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Arc.*mean$")))),
                                    rcp26_Arc_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Arc.*sd$")))^2)),
                                    rcp26_Ind_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Ind.*mean$")))),
                                    rcp26_Ind_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Ind.*sd$")))^2)),
                                    rcp26_Pac_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Pac.*mean$")))),
                                    rcp26_Pac_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp26_Pac.*sd$")))^2)),
                                    rcp85_Atl_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Atl.*mean$")))),
                                    rcp85_Atl_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Atl.*sd$")))^2)),
                                    rcp85_Arc_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Arc.*mean$")))),
                                    rcp85_Arc_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Arc.*sd$")))^2)),
                                    rcp85_Ind_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Ind.*mean$")))),
                                    rcp85_Ind_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Ind.*sd$")))^2)),
                                    rcp85_Pac_mean = rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Pac.*mean$")))),
                                    rcp85_Pac_sd = sqrt(rowMeans(as.data.frame(results_oceans %>% select(matches("^rcp85_Pac.*sd$")))^2)))

data.frame(control_Atl = mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp26_Atl_mean),mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp85_Atl_mean)),
           control_Atl_sd = sqrt(mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp26_Atl_sd^2),
                                      mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp85_Atl_sd^2))),
           control_Arc = mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp26_Arc_mean),mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp85_Arc_mean)),
           control_Arc_sd = sqrt(mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp26_Arc_sd^2),
                                      mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp85_Arc_sd^2))),
           control_Ind = mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp26_Ind_mean),mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp85_Ind_mean)),
           control_Ind_sd = sqrt(mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp26_Ind_sd^2),
                                      mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp85_Ind_sd^2))),
           control_Pac = mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp26_Pac_mean),mean(results_oceans_av_sd[results_oceans_av_sd$year <=2000,]$rcp85_Pac_mean)),
           control_Pac_sd = sqrt(mean(mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp26_Pac_sd^2),
                                      mean(results_oceans_av_sd[results_oceans_av_sd$year <= 2000,]$rcp85_Pac_sd^2))),
           
           future_Atl_rcp26 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Atl_mean),
           future_Atl_rcp26_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Atl_sd^2)),
           future_Arc_rcp26 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Arc_mean),
           future_Arc_rcp26_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Arc_sd^2)),
           future_Ind_rcp26 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Ind_mean),
           future_Ind_rcp26_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Ind_sd^2)),
           future_Pac_rcp26 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Pac_mean),
           future_Pac_rcp26_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp26_Pac_sd^2)),
           
           future_Atl_rcp85 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Atl_mean),
           future_Atl_rcp85_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Atl_sd^2)),
           future_Arc_rcp85 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Arc_mean),
           future_Arc_rcp85_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Arc_sd^2)),
           future_Ind_rcp85 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Ind_mean),
           future_Ind_rcp85_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Ind_sd^2)),
           future_Pac_rcp85 = mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Pac_mean),
           future_Pac_rcp85_sd = sqrt(mean(results_oceans_av_sd[results_oceans_av_sd$year >=2071,]$rcp85_Pac_sd^2))
)


# dataframe by climates
results_climats_av_sd = data.frame(year = results_climats$year,
                                 rcp26_A_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_A.*mean$")))),
                                 rcp26_A_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_A.*sd$")))^2)),
                                 rcp26_B_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_B.*mean$")))),
                                 rcp26_B_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_B.*sd$")))^2)),
                                 rcp26_C_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_C.*mean$")))),
                                 rcp26_C_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_C.*sd$")))^2)),
                                 rcp26_D_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_D.*mean$")))),
                                 rcp26_D_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_D.*sd$")))^2)),
                                 rcp26_E_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_E.*mean$")))),
                                 rcp26_E_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp26_E.*sd$")))^2)),
                                 rcp85_A_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_A.*mean$")))),
                                 rcp85_A_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_A.*sd$")))^2)),
                                 rcp85_B_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_B.*mean$")))),
                                 rcp85_B_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_B.*sd$")))^2)),
                                 rcp85_C_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_C.*mean$")))),
                                 rcp85_C_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_C.*sd$")))^2)),
                                 rcp85_D_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_D.*mean$")))),
                                 rcp85_D_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_D.*sd$")))^2)),
                                 rcp85_E_mean = rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_E.*mean$")))),
                                 rcp85_E_sd = sqrt(rowMeans(as.data.frame(results_climats %>% select(matches("^rcp85_E.*sd$")))^2)))

results_all = data.frame(year = results_climats_av_sd$year,
                         rcp26 = rowSums(as.data.frame(results_climats_av_sd %>% select(matches("^rcp26.*mean$")))),
                         rcp26_sd = rowSums(as.data.frame(results_climats_av_sd %>% select(matches("^rcp26.*sd$")))),
                         rcp85 = rowSums(as.data.frame(results_climats_av_sd %>% select(matches("^rcp85.*mean$")))),
                         rcp85_sd = rowSums(as.data.frame(results_climats_av_sd %>% select(matches("^rcp85.*sd$")))))


data.frame(control_global = mean(mean(results_all[results_all$year <=2000,]$rcp26),mean(results_all[results_all$year <=2000,]$rcp85)),
           control_global_sd = sqrt(mean(mean(results_all[results_all$year <= 2000,]$rcp26_sd^2),
                                         mean(results_all[results_all$year <= 2000,]$rcp85_sd^2))),
           future_rcp26 = mean(results_all[results_all$year >=2071,]$rcp26),
           future_rcp26_sd = sqrt(mean(results_all[results_all$year >=2071,]$rcp26_sd^2)),
           future_rcp85 = mean(results_all[results_all$year >=2071,]$rcp85),
           future_rcp85_sd = sqrt(mean(results_all[results_all$year >=2071,]$rcp85_sd^2)))

data.frame(control_A = mean(mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp26_A_mean),mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp85_A_mean)),
           control_A_sd = sqrt(mean(mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp26_A_sd^2),
                                    mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp85_A_sd^2))),
           control_B = mean(mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp26_B_mean),mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp85_B_mean)),
           control_B_sd = sqrt(mean(mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp26_B_sd^2),
                                    mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp85_B_sd^2))),
           control_C = mean(mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp26_C_mean),mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp85_C_mean)),
           control_C_sd = sqrt(mean(mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp26_C_sd^2),
                                    mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp85_C_sd^2))),
           control_D = mean(mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp26_D_mean),mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp85_D_mean)),
           control_D_sd = sqrt(mean(mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp26_D_sd^2),
                                    mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp85_D_sd^2))),
           control_E = mean(mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp26_E_mean),mean(results_climats_av_sd[results_climats_av_sd$year <=2000,]$rcp85_E_mean)),
           control_E_sd = sqrt(mean(mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp26_E_sd^2),
                                    mean(results_climats_av_sd[results_climats_av_sd$year <= 2000,]$rcp85_E_sd^2))),
           
           future_A_rcp26 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_A_mean),
           future_A_rcp26_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_A_sd^2)),
           future_B_rcp26 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_B_mean),
           future_B_rcp26_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_B_sd^2)),
           future_C_rcp26 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_C_mean),
           future_C_rcp26_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_C_sd^2)),
           future_D_rcp26 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_D_mean),
           future_D_rcp26_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_D_sd^2)),
           future_E_rcp26 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_E_mean),
           future_E_rcp26_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp26_E_sd^2)),
           
           future_A_rcp85 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_A_mean),
           future_A_rcp85_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_A_sd^2)),
           future_B_rcp85 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_B_mean),
           future_B_rcp85_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_B_sd^2)),
           future_C_rcp85 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_C_mean),
           future_C_rcp85_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_C_sd^2)),
           future_D_rcp85 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_D_mean),
           future_D_rcp85_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_D_sd^2)),
           future_E_rcp85 = mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_E_mean),
           future_E_rcp85_sd = sqrt(mean(results_climats_av_sd[results_climats_av_sd$year >=2071,]$rcp85_E_sd^2))
)


# Graphs ####
results_all$rcp26_color = ifelse(results_all$year <=2005, "black", "blue")
results_all$rcp85_color = ifelse(results_all$year <=2005, "black", "red")

results_climats$rcp26_color = ifelse(results_climats$year <=2005, "black", "blue")
results_climats$rcp85_color = ifelse(results_climats$year <=2005, "black", "red")

results_climats_av_sd$rcp26_color = ifelse(results_climats_av_sd$year <=2005, "black", "blue")
results_climats_av_sd$rcp85_color = ifelse(results_climats_av_sd$year <=2005, "black", "red")

results_oceans$rcp26_color = ifelse(results_oceans$year <=2005, "black", "blue")
results_oceans$rcp85_color = ifelse(results_oceans$year <=2005, "black", "red")

results_oceans_av_sd$rcp26_color = ifelse(results_oceans_av_sd$year <=2005, "black", "blue")
results_oceans_av_sd$rcp85_color = ifelse(results_oceans_av_sd$year <=2005, "black", "red")


list_rcp26_climats = c('rcp26_A', 'rcp26_B', 'rcp26_C', 'rcp26_D', 'rcp26_E')
list_rcp85_climats = c('rcp85_A', 'rcp85_B', 'rcp85_C', 'rcp85_D', 'rcp85_E')

list_rcp26_oceans = c('rcp26_Atl', 'rcp26_Arc', 'rcp26_Ind', 'rcp26_Pac')
list_rcp85_oceans = c('rcp85_Atl', 'rcp85_Arc', 'rcp85_Ind', 'rcp85_Pac')


# Trends
MC_trend = function(data, y, sd_y, N_MC){
  a=c(); b=c(); c=c(); mktau=c(); sensslope=c()
  for(i in 1:N_MC){
    pass_MC = data.frame(year = data$year, random_rcp= rnorm(nrow(data), mean=data[[y]], sd=data[[sd_y]]))
    a[i] = MannKendall(pass_MC$random_rcp)$sl
    mktau[i] = MannKendall(pass_MC$random_rcp)$tau
    b[i] = t.test(pass_MC[pass_MC$year <= 2000,]$random_rcp, pass_MC[pass_MC$year >= 2071,]$random_rcp)$p.value
    c[i] = sens.slope(pass_MC$random_rcp)$p.value
    sensslope[i] = sens.slope(pass_MC$random_rcp)$estimates[[1]]
  }
  return(data.frame(a=a, mktau=mktau, b=b, c=c, sensslope=sensslope))
}

rcp26_trend_MC = MC_trend(results_all, 'rcp26', 'rcp26_sd', 10000)
rcp85_trend_MC = MC_trend(results_all, 'rcp85', 'rcp85_sd', 10000)


summary_rcp_trends = format(data.frame(
  mktau= c(mean(rcp26_trend_MC$mktau),mean(rcp85_trend_MC$mktau)),
  mktau_std= c(sd(rcp26_trend_MC$mktau),sd(rcp85_trend_MC$mktau)),
  mannkendall=c(mean(rcp26_trend_MC$a), mean(rcp85_trend_MC$a)),
  mannkendall_std=c(sd(rcp26_trend_MC$a),sd(rcp85_trend_MC$a)),
  sensslope = c(mean(rcp26_trend_MC$sensslope),mean(rcp85_trend_MC$sensslope)),
  sensslope_std = c(sd(rcp26_trend_MC$sensslope),sd(rcp85_trend_MC$sensslope)),
  sen=c(mean(rcp26_trend_MC$c),mean(rcp85_trend_MC$c)),
  sen_std=c(sd(rcp26_trend_MC$c),sd(rcp85_trend_MC$c)),
  ttest=c(mean(rcp26_trend_MC$b),mean(rcp85_trend_MC$b)),
  ttest_std=c(sd(rcp26_trend_MC$b),sd(rcp85_trend_MC$b)),
  row.names = c('rcp26','rcp85')),digits=2)


rcp26_A_trend_MC = MC_trend(results_climats_av_sd, 'rcp26_A_mean', 'rcp26_A_sd', 10000)
rcp85_A_trend_MC = MC_trend(results_climats_av_sd, 'rcp85_A_mean', 'rcp85_A_sd', 10000)

rcp26_B_trend_MC = MC_trend(results_climats_av_sd, 'rcp26_B_mean', 'rcp26_B_sd', 10000)
rcp85_B_trend_MC = MC_trend(results_climats_av_sd, 'rcp85_B_mean', 'rcp85_B_sd', 10000)

rcp26_C_trend_MC = MC_trend(results_climats_av_sd, 'rcp26_C_mean', 'rcp26_C_sd', 10000)
rcp85_C_trend_MC = MC_trend(results_climats_av_sd, 'rcp85_C_mean', 'rcp85_C_sd', 10000)

rcp26_D_trend_MC = MC_trend(results_climats_av_sd, 'rcp26_D_mean', 'rcp26_D_sd', 10000)
rcp85_D_trend_MC = MC_trend(results_climats_av_sd, 'rcp85_D_mean', 'rcp85_D_sd', 10000)

rcp26_E_trend_MC = MC_trend(results_climats_av_sd, 'rcp26_E_mean', 'rcp26_E_sd', 10000)
rcp85_E_trend_MC = MC_trend(results_climats_av_sd, 'rcp85_E_mean', 'rcp85_E_sd', 10000)

rcp26_Atl_trend_MC = MC_trend(results_oceans_av_sd, 'rcp26_Atl_mean', 'rcp26_Atl_sd', 10000)
rcp85_Atl_trend_MC = MC_trend(results_oceans_av_sd, 'rcp85_Atl_mean', 'rcp85_Atl_sd', 10000)

rcp26_Arc_trend_MC = MC_trend(results_oceans_av_sd, 'rcp26_Arc_mean', 'rcp26_Arc_sd', 10000)
rcp85_Arc_trend_MC = MC_trend(results_oceans_av_sd, 'rcp85_Arc_mean', 'rcp85_Arc_sd', 10000)

rcp26_Ind_trend_MC = MC_trend(results_oceans_av_sd, 'rcp26_Ind_mean', 'rcp26_Ind_sd', 10000)
rcp85_Ind_trend_MC = MC_trend(results_oceans_av_sd, 'rcp85_Ind_mean', 'rcp85_Ind_sd', 10000)

rcp26_Pac_trend_MC = MC_trend(results_oceans_av_sd, 'rcp26_Pac_mean', 'rcp26_Pac_sd', 10000)
rcp85_Pac_trend_MC = MC_trend(results_oceans_av_sd, 'rcp85_Pac_mean', 'rcp85_Pac_sd', 10000)

summary_climats_trends = format(data.frame(mktau= c(mean(rcp26_A_trend_MC$mktau), mean(rcp26_B_trend_MC$mktau), mean(rcp26_C_trend_MC$mktau), mean(rcp26_D_trend_MC$mktau), mean(rcp26_E_trend_MC$mktau),
                                                      mean(rcp85_A_trend_MC$mktau), mean(rcp85_B_trend_MC$mktau), mean(rcp85_C_trend_MC$mktau), mean(rcp85_D_trend_MC$mktau), mean(rcp85_E_trend_MC$mktau)),
                                           mktau_std= c(sd(rcp26_A_trend_MC$mktau), sd(rcp26_B_trend_MC$mktau), sd(rcp26_C_trend_MC$mktau), sd(rcp26_D_trend_MC$mktau), sd(rcp26_E_trend_MC$mktau),
                                                          sd(rcp85_A_trend_MC$mktau), sd(rcp85_B_trend_MC$mktau), sd(rcp85_C_trend_MC$mktau), sd(rcp85_D_trend_MC$mktau), sd(rcp85_E_trend_MC$mktau)), 
                                           mannkendall=c(mean(rcp26_A_trend_MC$a), mean(rcp26_B_trend_MC$a), mean(rcp26_C_trend_MC$a), mean(rcp26_D_trend_MC$a), mean(rcp26_E_trend_MC$a),
                                                         mean(rcp85_A_trend_MC$a), mean(rcp85_B_trend_MC$a), mean(rcp85_C_trend_MC$a), mean(rcp85_D_trend_MC$a), mean(rcp85_E_trend_MC$a)),
                                           mannkendall_std=c(sd(rcp26_A_trend_MC$a), sd(rcp26_B_trend_MC$a), sd(rcp26_C_trend_MC$a), sd(rcp26_D_trend_MC$a), sd(rcp26_E_trend_MC$a),
                                                             sd(rcp85_A_trend_MC$a), sd(rcp85_B_trend_MC$a), sd(rcp85_C_trend_MC$a), sd(rcp85_D_trend_MC$a), sd(rcp85_E_trend_MC$a)),
                                           sensslope= c(mean(rcp26_A_trend_MC$sensslope), mean(rcp26_B_trend_MC$sensslope), mean(rcp26_C_trend_MC$sensslope), mean(rcp26_D_trend_MC$sensslope), mean(rcp26_E_trend_MC$sensslope),
                                                        mean(rcp85_A_trend_MC$sensslope), mean(rcp85_B_trend_MC$sensslope), mean(rcp85_C_trend_MC$sensslope), mean(rcp85_D_trend_MC$sensslope), mean(rcp85_E_trend_MC$sensslope)),
                                           sensslope_std= c(sd(rcp26_A_trend_MC$sensslope), sd(rcp26_B_trend_MC$sensslope), sd(rcp26_C_trend_MC$sensslope), sd(rcp26_D_trend_MC$sensslope), sd(rcp26_E_trend_MC$sensslope),
                                                            sd(rcp85_A_trend_MC$sensslope), sd(rcp85_B_trend_MC$sensslope), sd(rcp85_C_trend_MC$sensslope), sd(rcp85_D_trend_MC$sensslope), sd(rcp85_E_trend_MC$sensslope)), 
                                           sen=c(mean(rcp26_A_trend_MC$c), mean(rcp26_B_trend_MC$c), mean(rcp26_C_trend_MC$c), mean(rcp26_D_trend_MC$c), mean(rcp26_E_trend_MC$c),
                                                 mean(rcp85_A_trend_MC$c), mean(rcp85_B_trend_MC$c), mean(rcp85_C_trend_MC$c), mean(rcp85_D_trend_MC$c), mean(rcp85_E_trend_MC$c)),
                                           sen_std=c(sd(rcp26_A_trend_MC$c), sd(rcp26_B_trend_MC$c), sd(rcp26_C_trend_MC$c), sd(rcp26_D_trend_MC$c), sd(rcp26_E_trend_MC$c),
                                                     sd(rcp85_A_trend_MC$c), sd(rcp85_B_trend_MC$c), sd(rcp85_C_trend_MC$c), sd(rcp85_D_trend_MC$c), sd(rcp85_E_trend_MC$c)),
                                           ttest=c(mean(rcp26_A_trend_MC$b), mean(rcp26_B_trend_MC$b), mean(rcp26_C_trend_MC$b), mean(rcp26_D_trend_MC$b), mean(rcp26_E_trend_MC$b),
                                                   mean(rcp85_A_trend_MC$b), mean(rcp85_B_trend_MC$b), mean(rcp85_C_trend_MC$b), mean(rcp85_D_trend_MC$b), mean(rcp85_E_trend_MC$b)),
                                           ttest_std=c(sd(rcp26_A_trend_MC$b), sd(rcp26_B_trend_MC$b), sd(rcp26_C_trend_MC$b), sd(rcp26_D_trend_MC$b), sd(rcp26_E_trend_MC$b),
                                                       sd(rcp85_A_trend_MC$b), sd(rcp85_B_trend_MC$b), sd(rcp85_C_trend_MC$b), sd(rcp85_D_trend_MC$b), sd(rcp85_E_trend_MC$b)),
                                     row.names = c('rcp26_A','rcp26_B','rcp26_C','rcp26_D','rcp26_E',
                                                   'rcp85_A','rcp85_B','rcp85_C','rcp85_D','rcp85_E')), digits=2)

summary_oceans_trends = format(data.frame(mktau= c(mean(rcp26_Atl_trend_MC$mktau), mean(rcp26_Arc_trend_MC$mktau), mean(rcp26_Ind_trend_MC$mktau), mean(rcp26_Pac_trend_MC$mktau),
                                                      mean(rcp85_Atl_trend_MC$mktau), mean(rcp85_Arc_trend_MC$mktau), mean(rcp85_Ind_trend_MC$mktau), mean(rcp85_Pac_trend_MC$mktau)),
                                           mktau_std= c(sd(rcp26_Atl_trend_MC$mktau), sd(rcp26_Arc_trend_MC$mktau), sd(rcp26_Ind_trend_MC$mktau), sd(rcp26_Pac_trend_MC$mktau),
                                                          sd(rcp85_Atl_trend_MC$mktau), sd(rcp85_Arc_trend_MC$mktau), sd(rcp85_Ind_trend_MC$mktau), sd(rcp85_Pac_trend_MC$mktau)),
                                           mannkendall=c(mean(rcp26_Atl_trend_MC$a), mean(rcp26_Arc_trend_MC$a), mean(rcp26_Ind_trend_MC$a), mean(rcp26_Pac_trend_MC$a),
                                                         mean(rcp85_Atl_trend_MC$a), mean(rcp85_Arc_trend_MC$a), mean(rcp85_Ind_trend_MC$a), mean(rcp85_Pac_trend_MC$a)),
                                           mannkendall_std=c(sd(rcp26_Atl_trend_MC$a), sd(rcp26_Arc_trend_MC$a), sd(rcp26_Ind_trend_MC$a), sd(rcp26_Pac_trend_MC$a),
                                                             sd(rcp85_Atl_trend_MC$a), sd(rcp85_Arc_trend_MC$a), sd(rcp85_Ind_trend_MC$a), sd(rcp85_Pac_trend_MC$a)),
                                           sensslope= c(mean(rcp26_Atl_trend_MC$sensslope), mean(rcp26_Arc_trend_MC$sensslope), mean(rcp26_Ind_trend_MC$sensslope), mean(rcp26_Pac_trend_MC$sensslope),
                                                        mean(rcp85_Atl_trend_MC$sensslope), mean(rcp85_Arc_trend_MC$sensslope), mean(rcp85_Ind_trend_MC$sensslope), mean(rcp85_Pac_trend_MC$sensslope)),
                                           sensslope_std= c(sd(rcp26_Atl_trend_MC$sensslope), sd(rcp26_Arc_trend_MC$sensslope), sd(rcp26_Ind_trend_MC$sensslope), sd(rcp26_Pac_trend_MC$sensslope),
                                                            sd(rcp85_Atl_trend_MC$sensslope), sd(rcp85_Arc_trend_MC$sensslope), sd(rcp85_Ind_trend_MC$sensslope), sd(rcp85_Pac_trend_MC$sensslope)),
                                           sen=c(mean(rcp26_Atl_trend_MC$c), mean(rcp26_Arc_trend_MC$c), mean(rcp26_Ind_trend_MC$c), mean(rcp26_Pac_trend_MC$c),
                                                 mean(rcp85_Atl_trend_MC$c), mean(rcp85_Arc_trend_MC$c), mean(rcp85_Ind_trend_MC$c), mean(rcp85_Pac_trend_MC$c)),
                                           sen_std=c(sd(rcp26_Atl_trend_MC$c), sd(rcp26_Arc_trend_MC$c), sd(rcp26_Ind_trend_MC$c), sd(rcp26_Pac_trend_MC$c),
                                                     sd(rcp85_Atl_trend_MC$c), sd(rcp85_Arc_trend_MC$c), sd(rcp85_Ind_trend_MC$c), sd(rcp85_Pac_trend_MC$c)),
                                           ttest=c(mean(rcp26_Atl_trend_MC$b), mean(rcp26_Arc_trend_MC$b), mean(rcp26_Ind_trend_MC$b), mean(rcp26_Pac_trend_MC$b),
                                                   mean(rcp85_Atl_trend_MC$b), mean(rcp85_Arc_trend_MC$b), mean(rcp85_Ind_trend_MC$b), mean(rcp85_Pac_trend_MC$b)),
                                           ttest_std=c(sd(rcp26_Atl_trend_MC$b), sd(rcp26_Arc_trend_MC$b), sd(rcp26_Ind_trend_MC$b), sd(rcp26_Pac_trend_MC$b),
                                                       sd(rcp85_Atl_trend_MC$b), sd(rcp85_Arc_trend_MC$b), sd(rcp85_Ind_trend_MC$b), sd(rcp85_Pac_trend_MC$b)),
                                           row.names = c('rcp26_Atl','rcp26_Arc','rcp26_Ind','rcp26_Pac',
                                                         'rcp85_Atl','rcp85_Arc','rcp85_Ind','rcp85_Pac')), digits=2)

rcp26_trend_MC_Q = MC_trend(results_all_Q, 'rcp26', 'rcp26_sd', 10000)
rcp85_trend_MC_Q = MC_trend(results_all_Q, 'rcp85', 'rcp85_sd', 10000)

summary_rcp_trends_Q = format(data.frame(
  mktau= c(mean(rcp26_trend_MC_Q$mktau),mean(rcp85_trend_MC_Q$mktau)),
  mktau_std= c(sd(rcp26_trend_MC_Q$mktau),sd(rcp85_trend_MC_Q$mktau)),
  mannkendall=c(mean(rcp26_trend_MC_Q$a), mean(rcp85_trend_MC_Q$a)),
  mannkendall_std=c(sd(rcp26_trend_MC_Q$a),sd(rcp85_trend_MC_Q$a)),
  sensslope = c(mean(rcp26_trend_MC_Q$sensslope),mean(rcp85_trend_MC_Q$sensslope)),
  sensslope_std = c(sd(rcp26_trend_MC_Q$sensslope),sd(rcp85_trend_MC_Q$sensslope)),
  sen=c(mean(rcp26_trend_MC_Q$c),mean(rcp85_trend_MC_Q$c)),
  sen_std=c(sd(rcp26_trend_MC_Q$c),sd(rcp85_trend_MC_Q$c)),
  ttest=c(mean(rcp26_trend_MC_Q$b),mean(rcp85_trend_MC_Q$b)),
  ttest_std=c(sd(rcp26_trend_MC_Q$b),sd(rcp85_trend_MC_Q$b)),
  row.names = c('rcp26','rcp85')),digits=2)

rcp26_A_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp26_A_mean', 'rcp26_A_sd', 10000)
rcp85_A_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp85_A_mean', 'rcp85_A_sd', 10000)

rcp26_B_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp26_B_mean', 'rcp26_B_sd', 10000)
rcp85_B_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp85_B_mean', 'rcp85_B_sd', 10000)

rcp26_C_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp26_C_mean', 'rcp26_C_sd', 10000)
rcp85_C_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp85_C_mean', 'rcp85_C_sd', 10000)

rcp26_D_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp26_D_mean', 'rcp26_D_sd', 10000)
rcp85_D_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp85_D_mean', 'rcp85_D_sd', 10000)

rcp26_E_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp26_E_mean', 'rcp26_E_sd', 10000)
rcp85_E_trend_MC_Q = MC_trend(results_climats_Q_av_sd, 'rcp85_E_mean', 'rcp85_E_sd', 10000)

rcp26_Atl_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp26_Atl_mean', 'rcp26_Atl_sd', 10000)
rcp85_Atl_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp85_Atl_mean', 'rcp85_Atl_sd', 10000)

rcp26_Arc_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp26_Arc_mean', 'rcp26_Arc_sd', 10000)
rcp85_Arc_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp85_Arc_mean', 'rcp85_Arc_sd', 10000)

rcp26_Ind_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp26_Ind_mean', 'rcp26_Ind_sd', 10000)
rcp85_Ind_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp85_Ind_mean', 'rcp85_Ind_sd', 10000)

rcp26_Pac_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp26_Pac_mean', 'rcp26_Pac_sd', 10000)
rcp85_Pac_trend_MC_Q = MC_trend(results_oceans_Q_av_sd, 'rcp85_Pac_mean', 'rcp85_Pac_sd', 10000)

summary_climats_trends_Q = format(data.frame(mktau= c(mean(rcp26_A_trend_MC_Q$mktau), mean(rcp26_B_trend_MC_Q$mktau), mean(rcp26_C_trend_MC_Q$mktau), mean(rcp26_D_trend_MC_Q$mktau), mean(rcp26_E_trend_MC_Q$mktau),
                                                      mean(rcp85_A_trend_MC_Q$mktau), mean(rcp85_B_trend_MC_Q$mktau), mean(rcp85_C_trend_MC_Q$mktau), mean(rcp85_D_trend_MC_Q$mktau), mean(rcp85_E_trend_MC_Q$mktau)),
                                           mktau_std= c(sd(rcp26_A_trend_MC_Q$mktau), sd(rcp26_B_trend_MC_Q$mktau), sd(rcp26_C_trend_MC_Q$mktau), sd(rcp26_D_trend_MC_Q$mktau), sd(rcp26_E_trend_MC_Q$mktau),
                                                          sd(rcp85_A_trend_MC_Q$mktau), sd(rcp85_B_trend_MC_Q$mktau), sd(rcp85_C_trend_MC_Q$mktau), sd(rcp85_D_trend_MC_Q$mktau), sd(rcp85_E_trend_MC_Q$mktau)), 
                                           mannkendall=c(mean(rcp26_A_trend_MC_Q$a), mean(rcp26_B_trend_MC_Q$a), mean(rcp26_C_trend_MC_Q$a), mean(rcp26_D_trend_MC_Q$a), mean(rcp26_E_trend_MC_Q$a),
                                                         mean(rcp85_A_trend_MC_Q$a), mean(rcp85_B_trend_MC_Q$a), mean(rcp85_C_trend_MC_Q$a), mean(rcp85_D_trend_MC_Q$a), mean(rcp85_E_trend_MC_Q$a)),
                                           mannkendall_std=c(sd(rcp26_A_trend_MC_Q$a), sd(rcp26_B_trend_MC_Q$a), sd(rcp26_C_trend_MC_Q$a), sd(rcp26_D_trend_MC_Q$a), sd(rcp26_E_trend_MC_Q$a),
                                                             sd(rcp85_A_trend_MC_Q$a), sd(rcp85_B_trend_MC_Q$a), sd(rcp85_C_trend_MC_Q$a), sd(rcp85_D_trend_MC_Q$a), sd(rcp85_E_trend_MC_Q$a)),
                                           sensslope= c(mean(rcp26_A_trend_MC_Q$sensslope), mean(rcp26_B_trend_MC_Q$sensslope), mean(rcp26_C_trend_MC_Q$sensslope), mean(rcp26_D_trend_MC_Q$sensslope), mean(rcp26_E_trend_MC_Q$sensslope),
                                                        mean(rcp85_A_trend_MC_Q$sensslope), mean(rcp85_B_trend_MC_Q$sensslope), mean(rcp85_C_trend_MC_Q$sensslope), mean(rcp85_D_trend_MC_Q$sensslope), mean(rcp85_E_trend_MC_Q$sensslope)),
                                           sensslope_std= c(sd(rcp26_A_trend_MC_Q$sensslope), sd(rcp26_B_trend_MC_Q$sensslope), sd(rcp26_C_trend_MC_Q$sensslope), sd(rcp26_D_trend_MC_Q$sensslope), sd(rcp26_E_trend_MC_Q$sensslope),
                                                            sd(rcp85_A_trend_MC_Q$sensslope), sd(rcp85_B_trend_MC_Q$sensslope), sd(rcp85_C_trend_MC_Q$sensslope), sd(rcp85_D_trend_MC_Q$sensslope), sd(rcp85_E_trend_MC_Q$sensslope)), 
                                           sen=c(mean(rcp26_A_trend_MC_Q$c), mean(rcp26_B_trend_MC_Q$c), mean(rcp26_C_trend_MC_Q$c), mean(rcp26_D_trend_MC_Q$c), mean(rcp26_E_trend_MC_Q$c),
                                                 mean(rcp85_A_trend_MC_Q$c), mean(rcp85_B_trend_MC_Q$c), mean(rcp85_C_trend_MC_Q$c), mean(rcp85_D_trend_MC_Q$c), mean(rcp85_E_trend_MC_Q$c)),
                                           sen_std=c(sd(rcp26_A_trend_MC_Q$c), sd(rcp26_B_trend_MC_Q$c), sd(rcp26_C_trend_MC_Q$c), sd(rcp26_D_trend_MC_Q$c), sd(rcp26_E_trend_MC_Q$c),
                                                     sd(rcp85_A_trend_MC_Q$c), sd(rcp85_B_trend_MC_Q$c), sd(rcp85_C_trend_MC_Q$c), sd(rcp85_D_trend_MC_Q$c), sd(rcp85_E_trend_MC_Q$c)),
                                           ttest=c(mean(rcp26_A_trend_MC_Q$b), mean(rcp26_B_trend_MC_Q$b), mean(rcp26_C_trend_MC_Q$b), mean(rcp26_D_trend_MC_Q$b), mean(rcp26_E_trend_MC_Q$b),
                                                   mean(rcp85_A_trend_MC_Q$b), mean(rcp85_B_trend_MC_Q$b), mean(rcp85_C_trend_MC_Q$b), mean(rcp85_D_trend_MC_Q$b), mean(rcp85_E_trend_MC_Q$b)),
                                           ttest_std=c(sd(rcp26_A_trend_MC_Q$b), sd(rcp26_B_trend_MC_Q$b), sd(rcp26_C_trend_MC_Q$b), sd(rcp26_D_trend_MC_Q$b), sd(rcp26_E_trend_MC_Q$b),
                                                       sd(rcp85_A_trend_MC_Q$b), sd(rcp85_B_trend_MC_Q$b), sd(rcp85_C_trend_MC_Q$b), sd(rcp85_D_trend_MC_Q$b), sd(rcp85_E_trend_MC_Q$b)),
                                           row.names = c('rcp26_A','rcp26_B','rcp26_C','rcp26_D','rcp26_E',
                                                         'rcp85_A','rcp85_B','rcp85_C','rcp85_D','rcp85_E')), digits=2)

summary_oceans_trends_Q = format(data.frame(mktau= c(mean(rcp26_Atl_trend_MC_Q$mktau), mean(rcp26_Arc_trend_MC_Q$mktau), mean(rcp26_Ind_trend_MC_Q$mktau), mean(rcp26_Pac_trend_MC_Q$mktau),
                                                     mean(rcp85_Atl_trend_MC_Q$mktau), mean(rcp85_Arc_trend_MC_Q$mktau), mean(rcp85_Ind_trend_MC_Q$mktau), mean(rcp85_Pac_trend_MC_Q$mktau)),
                                          mktau_std= c(sd(rcp26_Atl_trend_MC_Q$mktau), sd(rcp26_Arc_trend_MC_Q$mktau), sd(rcp26_Ind_trend_MC_Q$mktau), sd(rcp26_Pac_trend_MC_Q$mktau),
                                                         sd(rcp85_Atl_trend_MC_Q$mktau), sd(rcp85_Arc_trend_MC_Q$mktau), sd(rcp85_Ind_trend_MC_Q$mktau), sd(rcp85_Pac_trend_MC_Q$mktau)),
                                          mannkendall=c(mean(rcp26_Atl_trend_MC_Q$a), mean(rcp26_Arc_trend_MC_Q$a), mean(rcp26_Ind_trend_MC_Q$a), mean(rcp26_Pac_trend_MC_Q$a),
                                                        mean(rcp85_Atl_trend_MC_Q$a), mean(rcp85_Arc_trend_MC_Q$a), mean(rcp85_Ind_trend_MC_Q$a), mean(rcp85_Pac_trend_MC_Q$a)),
                                          mannkendall_std=c(sd(rcp26_Atl_trend_MC_Q$a), sd(rcp26_Arc_trend_MC_Q$a), sd(rcp26_Ind_trend_MC_Q$a), sd(rcp26_Pac_trend_MC_Q$a),
                                                            sd(rcp85_Atl_trend_MC_Q$a), sd(rcp85_Arc_trend_MC_Q$a), sd(rcp85_Ind_trend_MC_Q$a), sd(rcp85_Pac_trend_MC_Q$a)),
                                          sensslope= c(mean(rcp26_Atl_trend_MC_Q$sensslope), mean(rcp26_Arc_trend_MC_Q$sensslope), mean(rcp26_Ind_trend_MC_Q$sensslope), mean(rcp26_Pac_trend_MC_Q$sensslope),
                                                       mean(rcp85_Atl_trend_MC_Q$sensslope), mean(rcp85_Arc_trend_MC_Q$sensslope), mean(rcp85_Ind_trend_MC_Q$sensslope), mean(rcp85_Pac_trend_MC_Q$sensslope)),
                                          sensslope_std= c(sd(rcp26_Atl_trend_MC_Q$sensslope), sd(rcp26_Arc_trend_MC_Q$sensslope), sd(rcp26_Ind_trend_MC_Q$sensslope), sd(rcp26_Pac_trend_MC_Q$sensslope),
                                                           sd(rcp85_Atl_trend_MC_Q$sensslope), sd(rcp85_Arc_trend_MC_Q$sensslope), sd(rcp85_Ind_trend_MC_Q$sensslope), sd(rcp85_Pac_trend_MC_Q$sensslope)),
                                          sen=c(mean(rcp26_Atl_trend_MC_Q$c), mean(rcp26_Arc_trend_MC_Q$c), mean(rcp26_Ind_trend_MC_Q$c), mean(rcp26_Pac_trend_MC_Q$c),
                                                mean(rcp85_Atl_trend_MC_Q$c), mean(rcp85_Arc_trend_MC_Q$c), mean(rcp85_Ind_trend_MC_Q$c), mean(rcp85_Pac_trend_MC_Q$c)),
                                          sen_std=c(sd(rcp26_Atl_trend_MC_Q$c), sd(rcp26_Arc_trend_MC_Q$c), sd(rcp26_Ind_trend_MC_Q$c), sd(rcp26_Pac_trend_MC_Q$c),
                                                    sd(rcp85_Atl_trend_MC_Q$c), sd(rcp85_Arc_trend_MC_Q$c), sd(rcp85_Ind_trend_MC_Q$c), sd(rcp85_Pac_trend_MC_Q$c)),
                                          ttest=c(mean(rcp26_Atl_trend_MC_Q$b), mean(rcp26_Arc_trend_MC_Q$b), mean(rcp26_Ind_trend_MC_Q$b), mean(rcp26_Pac_trend_MC_Q$b),
                                                  mean(rcp85_Atl_trend_MC_Q$b), mean(rcp85_Arc_trend_MC_Q$b), mean(rcp85_Ind_trend_MC_Q$b), mean(rcp85_Pac_trend_MC_Q$b)),
                                          ttest_std=c(sd(rcp26_Atl_trend_MC_Q$b), sd(rcp26_Arc_trend_MC_Q$b), sd(rcp26_Ind_trend_MC_Q$b), sd(rcp26_Pac_trend_MC_Q$b),
                                                      sd(rcp85_Atl_trend_MC_Q$b), sd(rcp85_Arc_trend_MC_Q$b), sd(rcp85_Ind_trend_MC_Q$b), sd(rcp85_Pac_trend_MC_Q$b)),
                                          row.names = c('rcp26_Atl','rcp26_Arc','rcp26_Ind','rcp26_Pac',
                                                        'rcp85_Atl','rcp85_Arc','rcp85_Ind','rcp85_Pac')), digits=2)


# Graph global ####
p1 = ggplot(data = results_all, aes(x = year)) +
  ggtitle('a') +
  xlab("") +
  ylab("Average DOC export (TgC "~y^-1*")") +
  expand_limits(y=0) +
  geom_line(aes(y=rcp26, color=rcp26_color), linewidth=1.2) +
  geom_ribbon(aes(x=year, ymin=rcp26-rcp26_sd, ymax=rcp26+rcp26_sd,
                  fill=rcp26_color), alpha=0.15, show.legend = FALSE) +
  geom_line(aes(y=rcp85, color=rcp85_color), linewidth=1.2) +
  geom_ribbon(aes(x=year, ymin=rcp85-rcp85_sd, ymax=rcp85+rcp85_sd,
                  fill=rcp85_color), alpha=0.15, show.legend = FALSE) +
  geom_rect(mapping=aes(xmin=1971, xmax=2000, ymin=85, ymax=170), colour = "black", alpha=0, linetype="dashed") +
  geom_rect(mapping=aes(xmin=2071, xmax=2099, ymin=90, ymax=200), colour = "black", alpha=0, linetype="dashed") +
  geom_text(mapping=aes(x=1977, y=175, label="Control"), color='black') +
  geom_text(mapping=aes(x=2077, y=206, label="Future"), color='black') +
  theme_publish() +
  scale_color_manual(values=c('black','blue','red'), name="RCP scenario",
                     labels=c("Historical","RCP 2.6","RCP 8.5")) +
  scale_fill_manual(values=c('black','blue','red'), name="RCP scenario",
                    labels=c("Historical","RCP 2.6","RCP 8.5")) +
  theme(legend.justification=c(0,1), legend.position=c(0.70,0.40),
        legend.background = element_rect(fill="white", linewidth=.5),
        legend.key.height=unit(2,"line"), legend.key.width=unit(2,"line"))


#graph diff
results_all$diff26 = results_all$rcp26-first(results_all$rcp26)
results_all$diff26_std = sqrt(abs(results_all$rcp26_sd^2-first(results_all$rcp26_sd^2)))
results_all$diff85 = results_all$rcp85-first(results_all$rcp85)
results_all$diff85_std = sqrt(abs(results_all$rcp85_sd^2-first(results_all$rcp85_sd^2)))

p2 = ggplot(data = results_all, aes(x = year)) +
  ggtitle('b') +
  xlab("") +
  ylab("Diff. in DOC export (TgC "~y^-1*")") +
  expand_limits(y=0) +
  geom_line(aes(y=diff26, color=rcp26_color), linewidth=1.2) +
  geom_ribbon(aes(x=year, ymin=diff26-diff26_std, ymax=diff26+diff26_std,
                  fill=rcp26_color), alpha=0.15, show.legend = FALSE) +
  geom_line(aes(y=diff85, color=rcp85_color), linewidth=1.2) +
  geom_ribbon(aes(x=year, ymin=diff85-diff85_std, ymax=diff85+diff85_std,
                  fill=rcp85_color), alpha=0.15, show.legend = FALSE) +
  geom_hline(yintercept=0) +
  geom_rect(mapping=aes(xmin=1971, xmax=2000, ymin=-21, ymax=23), colour = "black", alpha=0, linetype="dashed") +
  geom_rect(mapping=aes(xmin=2071, xmax=2099, ymin=-11, ymax=55), colour = "black", alpha=0, linetype="dashed") +
  theme_publish() +
  theme(legend.position="none") +
  scale_color_manual(values=c('black','blue','red'), name="RCP scenario",
                     labels=c("Historical","RCP 2.6","RCP 8.5")) +
  scale_fill_manual(values=c('black','blue','red'), name="RCP scenario",
                    labels=c("Historical","RCP 2.6","RCP 8.5"))

p26 = ggplot(rcp26_trend_MC) +
    xlab('p-value') +
    ylab("Frequency") +
    ggtitle('c') +
    geom_histogram(aes(a), binwidth = 0.01, fill='darkblue', alpha = 0.5) +
    geom_histogram(aes(b), binwidth = 0.01, fill='blue', alpha = 0.5) +
    geom_histogram(aes(c), binwidth = 0.01, fill='grey', alpha = 0.5) +
    geom_vline(xintercept=0.05, linetype="dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_publish()

p85 = ggplot(rcp85_trend_MC) +
  xlab('p-value') +
  ylab("Frequency") +
  ggtitle('d') +
  geom_histogram(aes(a), binwidth = 0.01, fill='darkred', alpha = 0.5) +
  geom_histogram(aes(b), binwidth = 0.01, fill='red', alpha = 0.5) +
  geom_histogram(aes(c), binwidth = 0.01, fill='grey', alpha = 0.5) +
  geom_vline(xintercept=0.05, linetype="dashed") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  theme_publish()

ggarrange(p1 | p2, p26 | p85, heights=c(0.8,0.2), ncol=1)






# Climate and Oceans Graphs ####
p1 = list()

for (i in 1:length(list_rcp26_climats)){
    pa = ggplot(data = results_climats_av_sd, aes(x = year)) +
      xlab("") +
      ylab("Average DOC export (TgC "~y^-1*")") +
      expand_limits(y=0) +
      ggtitle(all_climats_long[i]) +
      geom_line(aes(y=.data[[paste0(list_rcp26_climats[i],'_mean')]], color=rcp26_color), linewidth=1) +
      geom_ribbon(aes(x=year, ymin=.data[[paste0(list_rcp26_climats[i],'_mean')]]-.data[[paste0(list_rcp26_climats[i],'_sd')]], 
                      ymax=.data[[paste0(list_rcp26_climats[i],'_mean')]]+.data[[paste0(list_rcp26_climats[i],'_sd')]],
                      fill=rcp26_color), alpha=0.15, show.legend = FALSE) +
      geom_line(aes(y=.data[[paste0(list_rcp85_climats[i],'_mean')]], color=rcp85_color), linewidth=1) +
      geom_ribbon(aes(x=year, ymin=.data[[paste0(list_rcp85_climats[i],'_mean')]]-.data[[paste0(list_rcp85_climats[i],'_sd')]], 
                      ymax=.data[[paste0(list_rcp85_climats[i],'_mean')]]+.data[[paste0(list_rcp85_climats[i],'_sd')]],
                      fill=rcp85_color), alpha=0.15, show.legend = FALSE) +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
      theme_publish() +
      theme(legend.position="none") +
      scale_color_manual(values=c('black','blue','red'), name="RCP scenario",
                         labels=c("Historical","RCP 2.6","RCP 8.5")) +
      scale_fill_manual(values=c('black','blue','red'), name="RCP scenario",
                        labels=c("Historical","RCP 2.6","RCP 8.5"))
    pb = ggplot(get(paste0(list_rcp26_climats[i],'_trend_MC'))) +
      xlab('p-value') +
      ylab("Frequency") +
      expand_limits(x=c(0,1)) +
      geom_histogram(aes(a), binwidth = 0.01, fill='darkblue', alpha = 0.5) +
      geom_histogram(aes(b), binwidth = 0.01, fill='blue', alpha = 0.5) +
      geom_histogram(aes(c), binwidth = 0.01, fill='grey', alpha = 0.5) +
      geom_vline(xintercept=0.05, linetype="dashed") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
      theme_publish()
    
    pc = ggplot(get(paste0(list_rcp85_climats[i],'_trend_MC'))) +
      xlab('p-value') +
      ylab("Frequency") +
      expand_limits(x=c(0,1)) +
      geom_histogram(aes(a), binwidth = 0.01, fill='darkred', alpha = 0.5) +
      geom_histogram(aes(b), binwidth = 0.01, fill='red', alpha = 0.5) +
      geom_histogram(aes(c), binwidth = 0.01, fill='grey', alpha = 0.5) +
      geom_vline(xintercept=0.05, linetype="dashed") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
      theme_publish()
    p1[[i]] = list(pa, pb, pc)
}

p2 = list()

for (i in 1:length(list_rcp26_oceans)){
  pa = ggplot(data = results_oceans_av_sd, aes(x = year)) +
    xlab("") +
    ylab("Average DOC export (TgC "~y^-1*")") +
    expand_limits(y=0) +
    ggtitle(all_oceans[i]) +
    geom_line(aes(y=.data[[paste0(list_rcp26_oceans[i],'_mean')]], color=rcp26_color), linewidth=1) +
    geom_ribbon(aes(x=year, ymin=.data[[paste0(list_rcp26_oceans[i],'_mean')]]-.data[[paste0(list_rcp26_oceans[i],'_sd')]], 
                    ymax=.data[[paste0(list_rcp26_oceans[i],'_mean')]]+.data[[paste0(list_rcp26_oceans[i],'_sd')]],
                    fill=rcp26_color), alpha=0.15, show.legend = FALSE) +
    geom_line(aes(y=.data[[paste0(list_rcp85_oceans[i],'_mean')]], color=rcp85_color), linewidth=1) +
    geom_ribbon(aes(x=year, ymin=.data[[paste0(list_rcp85_oceans[i],'_mean')]]-.data[[paste0(list_rcp85_oceans[i],'_sd')]], 
                    ymax=.data[[paste0(list_rcp85_oceans[i],'_mean')]]+.data[[paste0(list_rcp85_oceans[i],'_sd')]],
                    fill=rcp85_color), alpha=0.15, show.legend = FALSE) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_publish() +
    theme(legend.position="none") +
    scale_color_manual(values=c('black','blue','red'), name="RCP scenario",
                       labels=c("Historical","RCP 2.6","RCP 8.5")) +
    scale_fill_manual(values=c('black','blue','red'), name="RCP scenario",
                      labels=c("Historical","RCP 2.6","RCP 8.5"))
  pb = ggplot(get(paste0(list_rcp26_oceans[i],'_trend_MC'))) +
    xlab('p-value') +
    ylab("Frequency") +
    expand_limits(x=c(0,1)) +
    geom_histogram(aes(a), binwidth = 0.01, fill='darkblue', alpha = 0.5) +
    geom_histogram(aes(b), binwidth = 0.01, fill='blue', alpha = 0.5) +
    geom_histogram(aes(c), binwidth = 0.01, fill='grey', alpha = 0.5) +
    geom_vline(xintercept=0.05, linetype="dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_publish()
  
  pc = ggplot(get(paste0(list_rcp85_oceans[i],'_trend_MC'))) +
    xlab('p-value') +
    ylab("Frequency") +
    expand_limits(x=c(0,1)) +
    geom_histogram(aes(a), binwidth = 0.01, fill='darkred', alpha = 0.5) +
    geom_histogram(aes(b), binwidth = 0.01, fill='red', alpha = 0.5) +
    geom_histogram(aes(c), binwidth = 0.01, fill='grey', alpha = 0.5) +
    geom_vline(xintercept=0.05, linetype="dashed") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
    theme_publish()
  p2[[i]] = list(pa, pb, pc)
}

ggarrange(
ggarrange((p1[[1]][[1]] / p1[[2]][[1]] / p1[[3]][[1]] / p1[[4]][[1]] / p1[[5]][[1]])+ plot_layout(axis_titles = "collect"),
          ((p1[[1]][[2]] / p1[[1]][[3]] / p1[[2]][[2]] / p1[[2]][[3]] / p1[[3]][[2]] / p1[[3]][[3]]/ 
              p1[[4]][[2]] / p1[[4]][[3]] / p1[[5]][[2]] / p1[[5]][[3]]) + plot_layout(axis_titles = "collect")), widths=c(1,1), ncol=2) |

ggarrange((p2[[2]][[1]] / p2[[1]][[1]] / p2[[3]][[1]] / p2[[4]][[1]] )+ plot_layout(axis_titles = "collect"),
          ((p2[[2]][[2]] / p2[[2]][[3]] / p2[[1]][[2]] / p2[[1]][[3]] / p2[[3]][[2]] / p2[[3]][[3]]/ 
              p2[[4]][[2]] / p2[[4]][[3]]) + plot_layout(axis_titles = "collect")), widths=c(1,1), ncol=2))





# daily graphs ####
# historical
daily_oceans_hist_av_sd = data.frame(julian = daily_oceans_hist$Julian,
                                rcp26_Atl_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Atl.*mean$")))),
                                rcp26_Atl_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Atl.*sd$")))^2)),
                                rcp26_Arc_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Arc.*mean$")))),
                                rcp26_Arc_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Arc.*sd$")))^2)),
                                rcp26_Ind_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Ind.*mean$")))),
                                rcp26_Ind_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Ind.*sd$")))^2)),
                                rcp26_Pac_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Pac.*mean$")))),
                                rcp26_Pac_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp26_Pac.*sd$")))^2)),
                                rcp85_Atl_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Atl.*mean$")))),
                                rcp85_Atl_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Atl.*sd$")))^2)),
                                rcp85_Arc_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Arc.*mean$")))),
                                rcp85_Arc_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Arc.*sd$")))^2)),
                                rcp85_Ind_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Ind.*mean$")))),
                                rcp85_Ind_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Ind.*sd$")))^2)),
                                rcp85_Pac_mean = rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Pac.*mean$")))),
                                rcp85_Pac_sd = sqrt(rowMeans(as.data.frame(daily_oceans_hist %>% select(matches("^rcp85_Pac.*sd$")))^2)))

daily_climats_hist_av_sd = data.frame(julian = daily_climats_hist$Julian,
                                 rcp26_A_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_A.*mean$")))),
                                 rcp26_A_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_A.*sd$")))^2)),
                                 rcp26_B_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_B.*mean$")))),
                                 rcp26_B_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_B.*sd$")))^2)),
                                 rcp26_C_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_C.*mean$")))),
                                 rcp26_C_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_C.*sd$")))^2)),
                                 rcp26_D_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_D.*mean$")))),
                                 rcp26_D_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_D.*sd$")))^2)),
                                 rcp26_E_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_E.*mean$")))),
                                 rcp26_E_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp26_E.*sd$")))^2)),
                                 rcp85_A_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_A.*mean$")))),
                                 rcp85_A_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_A.*sd$")))^2)),
                                 rcp85_B_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_B.*mean$")))),
                                 rcp85_B_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_B.*sd$")))^2)),
                                 rcp85_C_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_C.*mean$")))),
                                 rcp85_C_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_C.*sd$")))^2)),
                                 rcp85_D_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_D.*mean$")))),
                                 rcp85_D_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_D.*sd$")))^2)),
                                 rcp85_E_mean = rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_E.*mean$")))),
                                 rcp85_E_sd = sqrt(rowMeans(as.data.frame(daily_climats_hist %>% select(matches("^rcp85_E.*sd$")))^2)))

# future
daily_oceans_av_sd = data.frame(julian = daily_oceans$Julian,
                                  rcp26_Atl_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Atl.*mean$")))),
                                  rcp26_Atl_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Atl.*sd$")))^2)),
                                  rcp26_Arc_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Arc.*mean$")))),
                                  rcp26_Arc_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Arc.*sd$")))^2)),
                                  rcp26_Ind_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Ind.*mean$")))),
                                  rcp26_Ind_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Ind.*sd$")))^2)),
                                  rcp26_Pac_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Pac.*mean$")))),
                                  rcp26_Pac_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp26_Pac.*sd$")))^2)),
                                  rcp85_Atl_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Atl.*mean$")))),
                                  rcp85_Atl_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Atl.*sd$")))^2)),
                                  rcp85_Arc_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Arc.*mean$")))),
                                  rcp85_Arc_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Arc.*sd$")))^2)),
                                  rcp85_Ind_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Ind.*mean$")))),
                                  rcp85_Ind_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Ind.*sd$")))^2)),
                                  rcp85_Pac_mean = rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Pac.*mean$")))),
                                  rcp85_Pac_sd = sqrt(rowMeans(as.data.frame(daily_oceans %>% select(matches("^rcp85_Pac.*sd$")))^2)))

daily_climats_av_sd = data.frame(julian = daily_climats$Julian,
                                   rcp26_A_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_A.*mean$")))),
                                   rcp26_A_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_A.*sd$")))^2)),
                                   rcp26_B_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_B.*mean$")))),
                                   rcp26_B_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_B.*sd$")))^2)),
                                   rcp26_C_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_C.*mean$")))),
                                   rcp26_C_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_C.*sd$")))^2)),
                                   rcp26_D_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_D.*mean$")))),
                                   rcp26_D_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_D.*sd$")))^2)),
                                   rcp26_E_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_E.*mean$")))),
                                   rcp26_E_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp26_E.*sd$")))^2)),
                                   rcp85_A_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_A.*mean$")))),
                                   rcp85_A_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_A.*sd$")))^2)),
                                   rcp85_B_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_B.*mean$")))),
                                   rcp85_B_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_B.*sd$")))^2)),
                                   rcp85_C_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_C.*mean$")))),
                                   rcp85_C_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_C.*sd$")))^2)),
                                   rcp85_D_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_D.*mean$")))),
                                   rcp85_D_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_D.*sd$")))^2)),
                                   rcp85_E_mean = rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_E.*mean$")))),
                                   rcp85_E_sd = sqrt(rowMeans(as.data.frame(daily_climats %>% select(matches("^rcp85_E.*sd$")))^2)))


# DAILY CLIMATES


# daily_climats_diff_mean = 100*(daily_climats_av_sd[-1][ , c(TRUE,FALSE) ]-daily_climats_hist_av_sd[-1][ , c(TRUE,FALSE) ])/daily_climats_hist_av_sd[-1][ , c(TRUE,FALSE) ]
# 
# daily_climats_diff_sd = 100*sqrt((daily_climats_av_sd[-1][ , !c(TRUE,FALSE) ]^2*daily_climats_hist_av_sd[-1][ , c(TRUE,FALSE) ]^2+
#                                     daily_climats_hist_av_sd[-1][ , !c(TRUE,FALSE) ]^2*daily_climats_av_sd[-1][ , c(TRUE,FALSE) ]^2) /
#                                    daily_climats_hist_av_sd[-1][ , c(TRUE,FALSE) ]^4)

daily_climats_diff_mean = (daily_climats_av_sd[-1][ , c(TRUE,FALSE) ]-daily_climats_hist_av_sd[-1][ , c(TRUE,FALSE) ])
daily_climats_diff_sd = (daily_climats_av_sd[-1][ , !c(TRUE,FALSE) ]+daily_climats_hist_av_sd[-1][ , !c(TRUE,FALSE) ])

neworder <- order(c(2*(seq_along(daily_climats_diff_sd) - 1) + 1,
                    2*seq_along(daily_climats_diff_mean)))
daily_climats_diff = cbind(daily_climats_diff_mean, daily_climats_diff_sd)[,neworder]
daily_climats_diff$julian = daily_climats_hist_av_sd$julian

daily_oceans_diff_mean = (daily_oceans_av_sd[-1][ , c(TRUE,FALSE) ]-daily_oceans_hist_av_sd[-1][ , c(TRUE,FALSE) ])
daily_oceans_diff_sd = (daily_oceans_av_sd[-1][ , !c(TRUE,FALSE) ]+daily_oceans_hist_av_sd[-1][ , !c(TRUE,FALSE) ])

neworder <- order(c(2*(seq_along(daily_oceans_diff_sd) - 1) + 1,
                    2*seq_along(daily_oceans_diff_mean)))
daily_oceans_diff = cbind(daily_oceans_diff_mean, daily_oceans_diff_sd)[,neworder]
daily_oceans_diff$julian = daily_oceans_hist_av_sd$julian

p1_clim = ggplot(data = daily_climats_hist_av_sd, aes(x = julian)) +
  xlab("Julian day") +
  ylab("Av. DOC export (TgC."~d^-1*")") +
  ggtitle("Historical") +
  coord_cartesian(ylim=c(0,0.70)) +
  geom_line(aes(y=rcp26_A_mean, color="A"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_A_mean-rcp26_A_sd, ymax=rcp26_A_mean+rcp26_A_sd, fill="A"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_B_mean, color="B"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_B_mean-rcp26_B_sd, ymax=rcp26_B_mean+rcp26_B_sd, fill="B"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_C_mean, color="C"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_C_mean-rcp26_C_sd, ymax=rcp26_C_mean+rcp26_C_sd, fill="C"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_D_mean, color="D"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_D_mean-rcp26_D_sd, ymax=rcp26_D_mean+rcp26_D_sd, fill="D"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_E_mean, color="E"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_E_mean-rcp26_E_sd, ymax=rcp26_E_mean+rcp26_E_sd, fill="E"), alpha=0.3, show.legend = F) +
  theme_publish() +
  theme(legend.position = "None") +
  scale_color_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                    labels=c("Tropical","Semi-arid","Temperate","Cold","Polar")) +
  scale_fill_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                    labels=c("Tropical","Semi-arid","Temperate","Cold","Polar"))

p2_clim = ggplot(data = daily_climats_av_sd, aes(x = julian)) +
  xlab("") +
  ylab("") +
  ggtitle("RCP 2.6") +
  coord_cartesian(ylim=c(0,0.70)) +
  geom_line(aes(y=rcp26_A_mean, color="A"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_A_mean-rcp26_A_sd, ymax=rcp26_A_mean+rcp26_A_sd,
                  fill="A"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_B_mean, color="B"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_B_mean-rcp26_B_sd, ymax=rcp26_B_mean+rcp26_B_sd,
                  fill="B"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_C_mean, color="C"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_C_mean-rcp26_C_sd, ymax=rcp26_C_mean+rcp26_C_sd,
                  fill="C"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_D_mean, color="D"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_D_mean-rcp26_D_sd, ymax=rcp26_D_mean+rcp26_D_sd,
                  fill="D"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_E_mean, color="E"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_E_mean-rcp26_E_sd, ymax=rcp26_E_mean+rcp26_E_sd,
                  fill="E"), alpha=0.3, show.legend = F) +
  theme_publish() +
  theme(legend.position = "None") +
  scale_color_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                      labels=c("Tropical","Semi-arid","Temperate","Cold","Polar")) +
  scale_fill_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                    labels=c("Tropical","Semi-arid","Temperate","Cold","Polar"))

p3_clim = ggplot(data = daily_climats_av_sd, aes(x = julian)) +
  xlab("") +
  ylab("") +
  ggtitle("RCP 8.5") +
  coord_cartesian(ylim=c(0,0.70)) +
  geom_line(aes(y=rcp85_A_mean, color="A"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_A_mean-rcp85_A_sd, ymax=rcp85_A_mean+rcp85_A_sd,
                  fill="A"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_B_mean, color="B"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_B_mean-rcp85_B_sd, ymax=rcp85_B_mean+rcp85_B_sd,
                  fill="B"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_C_mean, color="C"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_C_mean-rcp85_C_sd, ymax=rcp85_C_mean+rcp85_C_sd,
                  fill="C"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_D_mean, color="D"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_D_mean-rcp85_D_sd, ymax=rcp85_D_mean+rcp85_D_sd,
                  fill="D"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_E_mean, color="E"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_E_mean-rcp85_E_sd, ymax=rcp85_E_mean+rcp85_E_sd,
                  fill="E"), alpha=0.3, show.legend = F) +
  theme_publish() +
  theme(legend.position = "None") +
  scale_color_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                      labels=c("Tropical","Semi-arid","Temperate","Cold","Polar")) +
  scale_fill_manual(values=c('darkgreen','goldenrod1','red','blue','turquoise'), name="Climate",
                    labels=c("Tropical","Semi-arid","Temperate","Cold","Polar"))

p4_clim = ggplot(data = daily_climats_diff, aes(x = julian)) +
  xlab("Julian day") +
  ylab(bquote(F[DOC]~change~(TgC~d^-1))) +
  coord_cartesian(ylim=c(-0.2,0.45)) +
  geom_line(aes(y=rcp26_A_mean, color="A"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_A_mean-rcp26_A_sd, ymax=rcp26_A_mean+rcp26_A_sd,
                  fill="A"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_B_mean, color="B"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_B_mean-rcp26_B_sd, ymax=rcp26_B_mean+rcp26_B_sd,
                  fill="B"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_C_mean, color="C"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_C_mean-rcp26_C_sd, ymax=rcp26_C_mean+rcp26_C_sd,
                  fill="C"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_D_mean, color="D"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_D_mean-rcp26_D_sd, ymax=rcp26_D_mean+rcp26_D_sd,
                  fill="D"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp26_E_mean, color="E"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_E_mean-rcp26_E_sd, ymax=rcp26_E_mean+rcp26_E_sd,
                  fill="E"), alpha=0.3, show.legend = F) +
  theme_publish() +
  scale_color_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                     labels=c("Tropical","Semi-arid","Temperate","Cold","Polar")) +
  scale_fill_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                    labels=c("Tropical","Semi-arid","Temperate","Cold","Polar"))

p5_clim = ggplot(data = daily_climats_diff, aes(x = julian)) +
  xlab("Julian day") +
  ylab("") +
  coord_cartesian(ylim=c(-0.2,0.45)) +
  geom_line(aes(y=rcp85_A_mean, color="A"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_A_mean-rcp85_A_sd, ymax=rcp85_A_mean+rcp85_A_sd,
                  fill="A"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_B_mean, color="B"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_B_mean-rcp85_B_sd, ymax=rcp85_B_mean+rcp85_B_sd,
                  fill="B"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_C_mean, color="C"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_C_mean-rcp85_C_sd, ymax=rcp85_C_mean+rcp85_C_sd,
                  fill="C"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_D_mean, color="D"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_D_mean-rcp85_D_sd, ymax=rcp85_D_mean+rcp85_D_sd,
                  fill="D"), alpha=0.3, show.legend = F) +
  geom_line(aes(y=rcp85_E_mean, color="E"), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_E_mean-rcp85_E_sd, ymax=rcp85_E_mean+rcp85_E_sd,
                  fill="E"), alpha=0.3, show.legend = F) +
  theme_publish() +
  theme(legend.position = "None") +
  scale_color_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                     labels=c("Tropical","Semi-arid","Temperate","Cold","Polar")) +
  scale_fill_manual(values=c("darkgreen","goldenrod1","red","blue","turquoise"), name="Climate",
                    labels=c("Tropical","Semi-arid","Temperate","Cold","Polar"))

# DAILY OCEANS
p1_oc = ggplot(data = daily_oceans_hist_av_sd, aes(x = julian)) +
  xlab("Julian day") +
  ylab("Av. DOC export (TgC."~d^-1*")") +
  ggtitle("Historical") +
  expand_limits(y=0) +
  geom_line(aes(y=rcp26_Arc_mean, color='Arc'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Arc_mean-rcp26_Arc_sd, ymax=rcp26_Arc_mean+rcp26_Arc_sd,
                  fill="Arc"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp26_Atl_mean, color='Atl'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Atl_mean-rcp26_Atl_sd, ymax=rcp26_Atl_mean+rcp26_Atl_sd,
                  fill="Atl"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp26_Ind_mean, color='Ind'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Ind_mean-rcp26_Ind_sd, ymax=rcp26_Ind_mean+rcp26_Ind_sd,
                  fill="Ind"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp26_Pac_mean, color='Pac'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Pac_mean-rcp26_Pac_sd, ymax=rcp26_Pac_mean+rcp26_Pac_sd,
                  fill="Pac"), alpha=0.3, show.legend = FALSE) +
  theme_publish() +
  theme(legend.position = "None") +
  coord_cartesian(ylim=c(0,0.60)) +
  scale_linetype_manual(name="Ocean supplied", values = c(rep("solid", 1), rep("solid", 1),
                                                          rep("dotted", 1), rep("dashed", 1)),
                        labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean")) +
  scale_color_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                     labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean")) +
  scale_fill_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                    labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean"))

p2_oc = ggplot(data = daily_oceans_av_sd, aes(x = julian)) +
  xlab("") +
  ylab("") +
  ggtitle("RCP 2.6") +
  expand_limits(y=0) +
  geom_line(aes(y=rcp26_Arc_mean, color='Arc'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Arc_mean-rcp26_Arc_sd, ymax=rcp26_Arc_mean+rcp26_Arc_sd,
                  fill="Arc"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp26_Atl_mean, color='Atl'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Atl_mean-rcp26_Atl_sd, ymax=rcp26_Atl_mean+rcp26_Atl_sd,
                  fill="Atl"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp26_Ind_mean, color='Ind'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Ind_mean-rcp26_Ind_sd, ymax=rcp26_Ind_mean+rcp26_Ind_sd,
                  fill="Ind"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp26_Pac_mean, color='Pac'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp26_Pac_mean-rcp26_Pac_sd, ymax=rcp26_Pac_mean+rcp26_Pac_sd,
                  fill="Pac"), alpha=0.3, show.legend = FALSE) +
  theme_publish() +
  theme(legend.position = "None") +
  coord_cartesian(ylim=c(0,0.60)) +
  scale_color_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                     labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean")) +
  scale_fill_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                    labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean"))

p3_oc = ggplot(data = daily_oceans_av_sd, aes(x = julian)) +
  xlab("") +
  ylab("") +
  ggtitle("RCP 8.5") +
  expand_limits(y=0) +
  geom_line(aes(y=rcp85_Arc_mean, color='Arc'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Arc_mean-rcp85_Arc_sd, ymax=rcp85_Arc_mean+rcp85_Arc_sd,
                  fill="Arc"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Atl_mean, color='Atl'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Atl_mean-rcp85_Atl_sd, ymax=rcp85_Atl_mean+rcp85_Atl_sd,
                  fill="Atl"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Ind_mean, color='Ind'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Ind_mean-rcp85_Ind_sd, ymax=rcp85_Ind_mean+rcp85_Ind_sd,
                  fill="Ind"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Pac_mean, color='Pac'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Pac_mean-rcp85_Pac_sd, ymax=rcp85_Pac_mean+rcp85_Pac_sd,
                  fill="Pac"), alpha=0.3, show.legend = FALSE) +
  theme_publish() +
  theme(legend.position = "None") +
  coord_cartesian(ylim=c(0,0.60)) +
  scale_color_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                     labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean")) +
  scale_fill_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                    labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean"))

p4_oc = ggplot(data = daily_oceans_diff, aes(x = julian)) +
  xlab("Julian day") +
  ylab(bquote(F[DOC]~change~(TgC~d^-1))) +
  coord_cartesian(ylim=c(-0.2,0.45)) +
  geom_line(aes(y=rcp85_Arc_mean, color='Arc'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Arc_mean-rcp85_Arc_sd, ymax=rcp85_Arc_mean+rcp85_Arc_sd,
                  fill="Arc"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Atl_mean, color='Atl'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Atl_mean-rcp85_Atl_sd, ymax=rcp85_Atl_mean+rcp85_Atl_sd,
                  fill="Atl"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Ind_mean, color='Ind'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Ind_mean-rcp85_Ind_sd, ymax=rcp85_Ind_mean+rcp85_Ind_sd,
                  fill="Ind"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Pac_mean, color='Pac'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Pac_mean-rcp85_Pac_sd, ymax=rcp85_Pac_mean+rcp85_Pac_sd,
                  fill="Pac"), alpha=0.3, show.legend = FALSE) +
  theme_publish() +
  scale_color_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                     labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean")) +
  scale_fill_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                    labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean"))

p5_oc = ggplot(data = daily_oceans_diff, aes(x = julian)) +
  xlab("Julian day") +
  ylab("") +
  coord_cartesian(ylim=c(-0.2,0.45)) +
  geom_line(aes(y=rcp85_Arc_mean, color='Arc'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Arc_mean-rcp85_Arc_sd, ymax=rcp85_Arc_mean+rcp85_Arc_sd,
                  fill="Arc"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Atl_mean, color='Atl'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Atl_mean-rcp85_Atl_sd, ymax=rcp85_Atl_mean+rcp85_Atl_sd,
                  fill="Atl"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Ind_mean, color='Ind'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Ind_mean-rcp85_Ind_sd, ymax=rcp85_Ind_mean+rcp85_Ind_sd,
                  fill="Ind"), alpha=0.3, show.legend = FALSE) +
  geom_line(aes(y=rcp85_Pac_mean, color='Pac'), linewidth=1.2) +
  geom_ribbon(aes(ymin=rcp85_Pac_mean-rcp85_Pac_sd, ymax=rcp85_Pac_mean+rcp85_Pac_sd,
                  fill="Pac"), alpha=0.3, show.legend = FALSE) +
  theme_publish() +
  theme(legend.position = "None") +
  scale_color_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                     labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean")) +
  scale_fill_manual(values=paletteer_d("colorBlindness::LightBlue2DarkBlue10Steps")[c(10,8,6,4)], name="Ocean supplied",
                    labels=c("Arctic Ocean","Atlantic Ocean","Indian Ocean","Pacific Ocean"))

ggarrange(ggarrange(((p1_clim/plot_spacer()|p2_clim/p4_clim|p3_clim/p5_clim) + plot_layout(axis_titles = "collect")), ncol=1),
          ggarrange(((p1_oc/plot_spacer()|p2_oc/p4_oc|p3_oc/p5_oc) + plot_layout(axis_titles = "collect")), ncol=1), ncol=1)

daily_oceans_hist_av_sd[daily_oceans_hist_av_sd$rcp26_Arc_mean==max(daily_oceans_hist_av_sd$rcp26_Arc_mean),]$julian
daily_oceans_hist_av_sd[daily_oceans_hist_av_sd$rcp85_Arc_mean==max(daily_oceans_hist_av_sd$rcp85_Arc_mean),]$julian
daily_oceans_av_sd[daily_oceans_av_sd$rcp26_Arc_mean==max(daily_oceans_av_sd$rcp26_Arc_mean),]$julian
daily_oceans_av_sd[daily_oceans_av_sd$rcp85_Arc_mean==max(daily_oceans_av_sd$rcp85_Arc_mean),]$julian
daily_oceans_diff[daily_oceans_diff$rcp26_Arc_mean==max(daily_oceans_diff$rcp26_Arc_mean),]$julian
daily_oceans_diff[daily_oceans_diff$rcp26_Arc_mean==min(daily_oceans_diff$rcp26_Arc_mean),]$julian


# dataframe by model
models_climats = format(data.frame(hist = colMeans(results_climats[1:(length(results_climats)-2)][results_climats[1:(length(results_climats)-2)]$year<=2000,]),
                          fut = colMeans(results_climats[1:(length(results_climats)-2)][results_climats[1:(length(results_climats)-2)]$year>=2071,])), digits=2)



models_oceans = format(data.frame(hist = colMeans(results_oceans[1:(length(results_oceans)-2)][results_oceans[1:(length(results_oceans)-2)]$year<=2000,]),
                                   fut = colMeans(results_oceans[1:(length(results_oceans)-2)][results_oceans[1:(length(results_oceans)-2)]$year>=2071,])), digits=2)





# Test if each catchment trend is ok ####
# path_catchments = paste0(dir,"/Final")
# 
# 
# outputs_p_values = expand.grid('DOC',RCPs,models,'p_value')
# outputs_p_values = outputs_p_values[order(outputs_p_values$Var1,outputs_p_values$Var2, outputs_p_values$Var3, outputs_p_values$Var4),]
# 
# climats[,paste(outputs_p_values$Var1,outputs_p_values$Var2,outputs_p_values$Var3,outputs_p_values$Var4, sep='_')] = NA
# 
# MC_trend_catch = function(data, y, sd_y, N_MC){
#   
#   cl = makeCluster(parallel::detectCores() - 1)
#   registerDoParallel(cl)
# 
# 
#   oper = foreach(l=1:N_MC, .packages=c('plyr','dplyr','Kendall','trend')) %dopar% {
#       pass_MC = data.frame(year = data$year, random_rcp= rnorm(nrow(data), mean=data[[y]], sd=data[[sd_y]]))
#       return(list(MannKendall(pass_MC$random_rcp)$sl, MannKendall(pass_MC$random_rcp)$tau, 
#                   t.test(pass_MC[pass_MC$year <= 2000,]$random_rcp, pass_MC[pass_MC$year >= 2071,]$random_rcp)$p.value, 
#                   sens.slope(pass_MC$random_rcp)$p.value, sens.slope(pass_MC$random_rcp)$estimates[[1]]))
#     }
# 
#     stopCluster(cl)
#     
#     result1 = do.call(rbind,lapply(oper,function(x){x[[1]]}))
#     result2 = do.call(rbind,lapply(oper,function(x){x[[2]]}))
#     result3 = do.call(rbind,lapply(oper,function(x){x[[3]]}))
#     result4 = do.call(rbind,lapply(oper,function(x){x[[4]]}))
#     result5 = do.call(rbind,lapply(oper,function(x){x[[5]]}))
#     
#     return(data.frame(a=result1, mktau=result2,
#                       b=result3, c=result4, sensslope=result5))
# }
# 
# for (name in row.names(climats)){
#   data_watershed = climats[row.names(climats)==name,]
#   print(c(name,data_watershed$CLIMAT,data_watershed$OCEAN))
#   area = data_watershed$Area * 1e6
# 
#   for (j in 1:length(RCPs)){
#     for (k in 1:length(models)){
#       print(c(models[k],RCPs[j]))
#       
#       file = read.csv(paste0(path_catchments,"/",name,"_",models[k],"_",RCPs[j],".csv"))
#       colnames(file)[colnames(file) == 'X'] <- 'year'
#       file_MC = MC_trend_catch(file, "mean", "sd", 10000)
# 
#       climats[row.names(climats)==name,][[paste0("DOC_",RCPs[j],"_",models[k],'_p_value')]] = mean(mean(file_MC$a),mean(file_MC$b),mean(file_MC$c))
#     }
#   }
# }
    
# write.csv(climats, file=file.path(path_final,paste0("climats_final.csv")))
  
#Maps ####
basins=read_sf(dsn = "E:/2016_2020_EcoLab/Couches/Global", layer = "basins_to_oceans")
basins <- st_zm(basins, drop=T)
basins <- st_make_valid(basins)
basins$diff_Q = NULL; basins$diff_Q_85 = NULL

climats$name_rivers = row.names(climats)

# if 3 out of 5 models significant
climats$p_value_26 = rowMedians(as.matrix(climats %>% select(matches("DOC_rcp26.*p_value$"))))
climats$p_value_85 = rowMedians(as.matrix(climats %>% select(matches("DOC_rcp85.*p_value$"))))


climats$tot_hist_26 = rowMeans(as.data.frame(climats %>% select(matches("DOC_hist_rcp26.*mean$"))))
climats$tot_hist_26_sd = rowMeans(as.data.frame(climats %>% select(matches("DOC_hist_rcp26.*sd$"))))
climats$tot_hist_85 = rowMeans(as.data.frame(climats %>% select(matches("DOC_hist_rcp85.*mean$"))))
climats$tot_hist_85_sd = rowMeans(as.data.frame(climats %>% select(matches("DOC_hist_rcp85.*mean$"))))
climats$tot_26 = rowMeans(as.data.frame(climats %>% select(matches("DOC_fut_rcp26.*mean$"))))
climats$tot_26_sd = rowMeans(as.data.frame(climats %>% select(matches("DOC_fut_rcp26.*sd$"))))
climats$tot_85 = rowMeans(as.data.frame(climats %>% select(matches("DOC_fut_rcp85.*mean$"))))
climats$tot_85_sd = rowMeans(as.data.frame(climats %>% select(matches("DOC_fut_rcp85.*mean$"))))

climats$diff_26 = ifelse(climats$p_value_26<0.05, 100* (climats$tot_26 - climats$tot_hist_26)/climats$tot_hist_26, 0)
climats$diff_85 = ifelse(climats$p_value_85<0.05, 100* (climats$tot_85 - climats$tot_hist_85)/climats$tot_hist_85, 0)

# climats$Q_hist_26 = rowMeans(data.frame(climats$Q_hist_rcp26_gfdl, climats$Q_hist_rcp26_hadgem, climats$Q_hist_rcp26_ipsl, climats$Q_hist_rcp26_miroc, climats$Q_hist_rcp26_noresm))
# climats$Q_hist_85 = rowMeans(data.frame(climats$Q_hist_rcp85_gfdl, climats$Q_hist_rcp85_hadgem, climats$Q_hist_rcp85_ipsl, climats$Q_hist_rcp85_miroc, climats$Q_hist_rcp85_noresm))
# climats$Q_26 = rowMeans(data.frame(climats$Q_rcp26_gfdl, climats$Q_rcp26_hadgem, climats$Q_rcp26_ipsl, climats$Q_rcp26_miroc, climats$Q_rcp26_noresm))
# climats$Q_85 = rowMeans(data.frame(climats$Q_rcp85_gfdl, climats$Q_rcp85_hadgem, climats$Q_rcp85_ipsl, climats$Q_rcp85_miroc, climats$Q_rcp85_noresm))
# climats$diff_Q_26 = 100* (climats$Q_26 - climats$Q_hist_26)/climats$Q_hist_26
# climats$diff_Q_85 = 100* (climats$Q_85 - climats$Q_hist_85)/climats$Q_hist_85
# 
# climats$SOC_hist_26 = rowMeans(data.frame(climats$SOC_hist_rcp26_gfdl, climats$SOC_hist_rcp26_hadgem, climats$SOC_hist_rcp26_ipsl, climats$SOC_hist_rcp26_miroc, climats$SOC_hist_rcp26_noresm))
# climats$SOC_hist_85 = rowMeans(data.frame(climats$SOC_hist_rcp85_gfdl, climats$SOC_hist_rcp85_hadgem, climats$SOC_hist_rcp85_ipsl, climats$SOC_hist_rcp85_miroc, climats$SOC_hist_rcp85_noresm))
# climats$SOC_26 = rowMeans(data.frame(climats$SOC_rcp26_gfdl, climats$SOC_rcp26_hadgem, climats$SOC_rcp26_ipsl, climats$SOC_rcp26_miroc, climats$SOC_rcp26_noresm))
# climats$SOC_85 = rowMeans(data.frame(climats$SOC_rcp85_gfdl, climats$SOC_rcp85_hadgem, climats$SOC_rcp85_ipsl, climats$SOC_rcp85_miroc, climats$SOC_rcp85_noresm))
# climats$diff_SOC_26 = 100* (climats$SOC_26 - climats$SOC_hist_26)/climats$SOC_hist_26
# climats$diff_SOC_85 = 100* (climats$SOC_85 - climats$SOC_hist_85)/climats$SOC_hist_85
# 
# climats$PCP_hist_26 = rowMeans(data.frame(climats$PCP_hist_rcp26_gfdl, climats$PCP_hist_rcp26_hadgem, climats$PCP_hist_rcp26_ipsl, climats$PCP_hist_rcp26_miroc, climats$PCP_hist_rcp26_noresm))
# climats$PCP_hist_85 = rowMeans(data.frame(climats$PCP_hist_rcp85_gfdl, climats$PCP_hist_rcp85_hadgem, climats$PCP_hist_rcp85_ipsl, climats$PCP_hist_rcp85_miroc, climats$PCP_hist_rcp85_noresm))
# climats$PCP_26 = rowMeans(data.frame(climats$PCP_rcp26_gfdl, climats$PCP_rcp26_hadgem, climats$PCP_rcp26_ipsl, climats$PCP_rcp26_miroc, climats$PCP_rcp26_noresm))
# climats$PCP_85 = rowMeans(data.frame(climats$PCP_rcp85_gfdl, climats$PCP_rcp85_hadgem, climats$PCP_rcp85_ipsl, climats$PCP_rcp85_miroc, climats$PCP_rcp85_noresm))
# climats$diff_PCP_26 = 100* (climats$PCP_26 - climats$PCP_hist_26)/climats$PCP_hist_26
# climats$diff_PCP_85 = 100* (climats$PCP_85 - climats$PCP_hist_85)/climats$PCP_hist_85
# climats$dev_PCP_26 = climats$PCP_26 - climats$PCP_hist_26
# climats$dev_PCP_85 = climats$PCP_85 - climats$PCP_hist_85
# 
# climats$TMP_hist_26 = rowMeans(data.frame(climats$TMP_hist_rcp26_gfdl, climats$TMP_hist_rcp26_hadgem, climats$TMP_hist_rcp26_ipsl, climats$TMP_hist_rcp26_miroc, climats$TMP_hist_rcp26_noresm))
# climats$TMP_hist_85 = rowMeans(data.frame(climats$TMP_hist_rcp85_gfdl, climats$TMP_hist_rcp85_hadgem, climats$TMP_hist_rcp85_ipsl, climats$TMP_hist_rcp85_miroc, climats$TMP_hist_rcp85_noresm))
# climats$TMP_26 = rowMeans(data.frame(climats$TMP_rcp26_gfdl, climats$TMP_rcp26_hadgem, climats$TMP_rcp26_ipsl, climats$TMP_rcp26_miroc, climats$TMP_rcp26_noresm))
# climats$TMP_85 = rowMeans(data.frame(climats$TMP_rcp85_gfdl, climats$TMP_rcp85_hadgem, climats$TMP_rcp85_ipsl, climats$TMP_rcp85_miroc, climats$TMP_rcp85_noresm))
# climats$diff_TMP_26 = 100* (climats$TMP_26 - climats$TMP_hist_26)/abs(climats$TMP_hist_26)
# climats$diff_TMP_85 = 100* (climats$TMP_85 - climats$TMP_hist_85)/abs(climats$TMP_hist_85)
# climats$dev_TMP_26 = climats$TMP_26 - climats$TMP_hist_26
# climats$dev_TMP_85 = climats$TMP_85 - climats$TMP_hist_85

basins_map=merge(basins, climats, by.x = "wribasin_u", by.y = "name_rivers")

basins_map$diff_26_cut = cut(basins_map$diff_26, breaks=c(-75,-50,-10,-5,-0.1,0.1,5,10,50,100,650), right = FALSE)
basins_map$diff_85_cut = cut(basins_map$diff_85, breaks=c(-75,-50,-10,-5,-0.1,0.1,5,10,50,100,650), right = FALSE)

palette_color = c('#b2182b','#d6604d','#f4a582','#fddbc7','lightgrey','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')
data("World")
pmap26 = ggplot(World) +
          geom_sf(color='transparent') +
          geom_sf(data=basins_map, aes(fill = diff_26_cut)
                  , color = "black", linewidth = .4) +
          scale_fill_manual(name="DOC change (%)",
                            labels=c("< -50", "-50 to -10", "-10 to -5", "-5 to 0","0", "0 to 5", "5 to 10", "10 to 50", "50 to 100", "> 100"),
                            values = palette_color,
                            drop = FALSE) +
          theme_publish() +
          theme(plot.margin = unit(c(0, 0, 0, 0), "null"))

pmap85 = ggplot(World) +
  geom_sf(color='transparent') +
  geom_sf(data=basins_map, aes(fill = diff_85_cut)
          , color = "black", linewidth = .4) +
  scale_fill_manual(name="DOC change (%)",
                    labels=c("< -50", "-50 to -10", "-10 to -5", "-5 to 0","0", "0 to 5", "5 to 10", "10 to 50", "50 to 100", "> 100"),
                    values = palette_color,
                    drop = FALSE) +
  theme_publish() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "null"))

# g_legend<-function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
# 
# mylegend<-g_legend(pmap85)



climats$group_lat = cut(climats$LAT_mouth, breaks=seq(-90,90,10), right = FALSE)
climats$group_lon = cut(climats$LONG_mouth, breaks=seq(-180,180,20), right = FALSE)

DOC_lat = as.data.frame(ddply(climats,~group_lat,summarise,
                              DOC_26_hist = mean(sum(DOC_hist_rcp26_gfdl_mean),
                                                 sum(DOC_hist_rcp26_hadgem_mean),
                                                 sum(DOC_hist_rcp26_ipsl_mean),
                                                 sum(DOC_hist_rcp26_miroc_mean),
                                                 sum(DOC_hist_rcp26_noresm_mean)),
                              DOC_26_hist_sd = sqrt(mean(sum(DOC_hist_rcp26_gfdl_sd^2),
                                                    sum(DOC_hist_rcp26_hadgem_sd^2),
                                                    sum(DOC_hist_rcp26_ipsl_sd^2),
                                                    sum(DOC_hist_rcp26_miroc_sd^2),
                                                    sum(DOC_hist_rcp26_noresm_sd^2))),
                              DOC_26_fut = mean(sum(DOC_fut_rcp26_gfdl_mean),
                                                sum(DOC_fut_rcp26_hadgem_mean),
                                                sum(DOC_fut_rcp26_ipsl_mean),
                                                sum(DOC_fut_rcp26_miroc_mean),
                                                sum(DOC_fut_rcp26_noresm_mean)),
                              DOC_26_fut_sd = sqrt(mean(sum(DOC_fut_rcp26_gfdl_sd^2),
                                                         sum(DOC_fut_rcp26_hadgem_sd^2),
                                                         sum(DOC_fut_rcp26_ipsl_sd^2),
                                                         sum(DOC_fut_rcp26_miroc_sd^2),
                                                         sum(DOC_fut_rcp26_noresm_sd^2))),
                              DOC_85_hist = mean(sum(DOC_hist_rcp85_gfdl_mean),
                                                 sum(DOC_hist_rcp85_hadgem_mean),
                                                 sum(DOC_hist_rcp85_ipsl_mean),
                                                 sum(DOC_hist_rcp85_miroc_mean),
                                                 sum(DOC_hist_rcp85_noresm_mean)),
                              DOC_85_hist_sd = sqrt(mean(sum(DOC_hist_rcp85_gfdl_sd^2),
                                                         sum(DOC_hist_rcp85_hadgem_sd^2),
                                                         sum(DOC_hist_rcp85_ipsl_sd^2),
                                                         sum(DOC_hist_rcp85_miroc_sd^2),
                                                         sum(DOC_hist_rcp85_noresm_sd^2))),
                              DOC_85_fut = mean(sum(DOC_fut_rcp85_gfdl_mean),
                                                sum(DOC_fut_rcp85_hadgem_mean),
                                                sum(DOC_fut_rcp85_ipsl_mean),
                                                sum(DOC_fut_rcp85_miroc_mean),
                                                sum(DOC_fut_rcp85_noresm_mean)),
                              DOC_85_fut_sd = sqrt(mean(sum(DOC_fut_rcp85_gfdl_sd^2),
                                                         sum(DOC_fut_rcp85_hadgem_sd^2),
                                                         sum(DOC_fut_rcp85_ipsl_sd^2),
                                                         sum(DOC_fut_rcp85_miroc_sd^2),
                                                         sum(DOC_fut_rcp85_noresm_sd^2)))))



DOC_lat$p_value_26 = NA; DOC_lat$p_value_85 = NA

for(i in 1:nrow(DOC_lat)){
    test_26 = data.frame(random_rcp_hist= rnorm(10000, DOC_lat$DOC_26_hist[i], DOC_lat$DOC_26_hist_sd[i]),
                         random_rcp_fut= rnorm(10000, DOC_lat$DOC_26_fut[i], DOC_lat$DOC_26_fut_sd[i]))
    test_85 = data.frame(random_rcp_hist= rnorm(10000, DOC_lat$DOC_85_hist[i], DOC_lat$DOC_85_hist_sd[i]),
                         random_rcp_fut= rnorm(10000, DOC_lat$DOC_85_fut[i], DOC_lat$DOC_85_fut_sd[i])) 
print(t.test(test_26$random_rcp_hist, test_26$random_rcp_fut)$p.value)
    DOC_lat[i,]['p_value_26'] = t.test(test_26$random_rcp_hist, test_26$random_rcp_fut)$p.value
    DOC_lat[i,]['p_value_85'] = t.test(test_85$random_rcp_hist, test_85$random_rcp_fut)$p.value
}


DOC_lat$change_26 = 100*(DOC_lat$DOC_26_fut - DOC_lat$DOC_26_hist) / DOC_lat$DOC_26_hist
# 100 * sqrt((var(y)*mean(x)^2 + var(x)*mean(y)^2)/(mean(y)^4)
DOC_lat$change_26_sd = 100*sqrt((DOC_lat$DOC_26_fut_sd^2*DOC_lat$DOC_26_hist^2+
                                  DOC_lat$DOC_26_hist_sd^2*DOC_lat$DOC_26_fut^2)/
                                  DOC_lat$DOC_26_hist^4)

DOC_lat$change_85 = 100*(DOC_lat$DOC_85_fut - DOC_lat$DOC_85_hist) / DOC_lat$DOC_85_hist
DOC_lat$change_85_sd = 100*sqrt((DOC_lat$DOC_85_fut_sd^2*DOC_lat$DOC_85_hist^2+
                                   DOC_lat$DOC_85_hist_sd^2*DOC_lat$DOC_85_fut^2)/
                                  DOC_lat$DOC_85_hist^4)

# DOC_lon = as.data.frame(ddply(climats,~group_lon,summarise,
#                               DOC_26_hist = mean(sum(DOC_hist_rcp26_gfdl_mean),
#                                                  sum(DOC_hist_rcp26_hadgem_mean),
#                                                  sum(DOC_hist_rcp26_ipsl_mean),
#                                                  sum(DOC_hist_rcp26_miroc_mean),
#                                                  sum(DOC_hist_rcp26_noresm_mean)),
#                               DOC_26_fut = mean(sum(DOC_fut_rcp26_gfdl_mean),
#                                                 sum(DOC_fut_rcp26_hadgem_mean),
#                                                 sum(DOC_fut_rcp26_ipsl_mean),
#                                                 sum(DOC_fut_rcp26_miroc_mean),
#                                                 sum(DOC_fut_rcp26_noresm_mean)),
#                               DOC_85_hist = mean(sum(DOC_hist_rcp85_gfdl_mean),
#                                                  sum(DOC_hist_rcp85_hadgem_mean),
#                                                  sum(DOC_hist_rcp85_ipsl_mean),
#                                                  sum(DOC_hist_rcp85_miroc_mean),
#                                                  sum(DOC_hist_rcp85_noresm_mean)),
#                               DOC_85_fut = mean(sum(DOC_fut_rcp85_gfdl_mean),
#                                                 sum(DOC_fut_rcp85_hadgem_mean),
#                                                 sum(DOC_fut_rcp85_ipsl_mean),
#                                                 sum(DOC_fut_rcp85_miroc_mean),
#                                                 sum(DOC_fut_rcp85_noresm_mean))))
# 
# DOC_lon$change_26 = 100*(DOC_lon$DOC_26_fut - DOC_lon$DOC_26_hist) / DOC_lon$DOC_26_hist
# DOC_lon$change_85 = 100*(DOC_lon$DOC_85_fut - DOC_lon$DOC_85_hist) / DOC_lon$DOC_85_hist

latitude=c('50-60S','40-50S','30-40S','20-30S','10-20S','0-10S',
           '0-10N','10-20N','20-30N','30-40N','40-50N','50-60N','60-70N','70-80N')
longitude=c('50-60S','40-50S','30-40S','20-30S','10-20S','0-10S',
            '0-10N','10-20N','20-30N','30-40N','40-50N','50-60N','60-70N','70-80N')

plat26 = ggplot(DOC_lat, aes(x=group_lat, y=change_26)) +
  geom_bar(stat="identity", color="black", fill="white", position=position_dodge()) +
  geom_errorbar(aes(ymin=change_26-change_26_sd, ymax=change_26+change_26_sd), width=.2,
                position=position_dodge(.9)) +
  theme_publish() +
  ylim(-65,110) +
  scale_x_discrete(labels= latitude) +
  coord_flip() +
  labs(x="", y = "")

plat85 = ggplot(DOC_lat, aes(x=group_lat, y=change_85)) +
  geom_bar(stat="identity", color="black",fill="white", position=position_dodge()) +
  geom_errorbar(aes(ymin=change_85-change_85_sd, ymax=change_85+change_85_sd), width=.2,
                position=position_dodge(.9)) +
  theme_publish() +
  ylim(-65,110) +
  scale_x_discrete(labels= latitude) +
  coord_flip() +
  labs(x="", y = "")

ggarrange(pmap26, plat26, pmap85, plat85, ncol = 2, nrow = 2, widths=c(2,1))

# ### Map total
# map_26_past = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="tot_hist_26", breaks=c(0,0.5,1.5,3.5,16,25),
#           title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
#           palette="YlOrRd") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_26_past, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/map_26_past.png")
# 
# map_85_past = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="tot_hist_85", breaks=c(0,0.5,1.5,3.5,16,25),
#           title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
#           palette="YlOrRd") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_85_past, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/map_85_past.png")
# 
# map_26_fut = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="tot_26", breaks=c(0,0.5,1.5,3.5,16,25),
#           title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
#           palette="YlOrRd") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_26_fut, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/map_26_fut.png")
# 
# map_85_fut = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="tot_85", breaks=c(0,0.5,1.5,3.5,16,25),
#           title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
#           palette="YlOrRd") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_85_fut, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/map_85_fut.png")


# Map Q, ORGC, Precip & Temp ####

# Discharge ####
# map_Q_hist = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="Q_hist_26", breaks=c(0,100,200,500,1000,10000),
#           labels=c("< 100", "100 to 200", "200 to 500", "500 to 1000", "> 1000"),
#           title=expression("Discharge ("~km^3*""~yr^-1*")"),
#           palette="Blues") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_Q_hist, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/tot_Q_map_past.png")
# 
# map_Q_26 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="diff_Q_26", breaks=c(-100,-10,-5,0,5,10,1000),
#           labels=c("< -10","-10 to -5","-5 to 0",
#                    "0 to 5","5 to 10","> 10"),
#           title=expression("Discharge change (%)"),
#           palette="RdBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_Q_26, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/change_Q_map_26.png")
# 
# map_Q_85 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="diff_Q_85", breaks=c(-100,-10,-5,0,5,10,1000),
#           labels=c("< -10","-10 to -5","-5 to 0",
#                    "0 to 5","5 to 10","> 10"),
#           title=expression("Discharge change (%)"),
#           palette="RdBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_Q_85, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/change_Q_map_85.png")
# 
# 
# # Soil C ####
# map_soil_hist = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="SOC_hist_26", breaks=c(0,1,5,10,25,50,100),
#           labels=c("< 1", "1 to 5","5 to 10", "10 to 25", "25 to 50", "> 50"),
#           title=expression("Soil organic carbon (g"~m^-3*")"),
#           palette="YlOrBr") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_soil_hist, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/soil_map_past.png")
# 
# map_soil_26 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="diff_SOC_26", breaks=c(-100,-50,-10,0,10,50,100),
#           labels=c("< -50","-50 to -10","-10 to 0",
#                    "0 to 10","10 to 50","> 50"),
#           title=expression("SOC change (%)"),
#           palette="RdYlBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_soil_26, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/soil_map_26.png")
# 
# map_soil_85 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="diff_SOC_85", breaks=c(-100,-50,-10,0,10,50,100),
#           labels=c("< -50","-50 to -10","-10 to 0",
#                    "0 to 10","10 to 50","> 50"),
#           title=expression("SOC change (%)"),
#           palette="RdYlBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_soil_85, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/soil_map_85.png")
# 
# 
# # Temperature ####
# map_T_past = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="TMP_hist_26", breaks=c(-20,-10,-5,0,5,10,15,20,25,100),
#           labels=c("< -10", "-10 to -5", "-5 to 0", "0 to 5", "5 to 10", "10 to 15",
#                    "15 to 20", "20 to 25", "> 25"),
#           title=expression("Air temperature (°C)"),
#           palette="-RdYlBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_T_past, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/T_map_past.png")
# 
# map_T_26 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="dev_TMP_26", breaks=c(-15,0,1,2,5,10),
#           labels=c("-15 to 0", "0 to 1", "1 to 2", "2 to 5", "5 to 10"),
#           title=expression("Temperature deviation ("~Delta*"°C)"),
#           palette="-RdYlBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_T_26, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/T_map_26.png")
# 
# map_T_85 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="dev_TMP_85", breaks=c(-15,0,1,2,5,10),
#           labels=c("-15 to 0", "0 to 1", "1 to 2", "2 to 5", "5 to 10"),
#           title=expression("Temperature deviation ("~Delta*"°C)"),
#           palette="-RdYlBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_T_85, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/T_map_85.png")
# 
# # Precipitation ####
# map_P_past = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="PCP_hist_26", breaks=c(0,200,400,600,800,1000,10000),
#           labels=c("0 to 200", "200 to 400", "400 to 600", "600 to 800", "800 to 1000", "> 1000"),
#           title=expression("Precipitation (mm"~yr^-1*")"),
#           palette="Blues") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.5,
#             legend.text.size = 1.2)
# tmap_save(tm = map_P_past, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/P_map_past.png")
# 
# map_P_26 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="dev_PCP_26", breaks=c(-600,-200,-100,0,100,200,500,1000),
#           labels=c("< -200", "-200 to -100", "-100 to 0", "0 to 100", "100 to 200", "200 to 500", "> 500"),
#           title=expression("Precipitation deviation (mm"~yr^-1*")"),
#           palette="RdBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.3,
#             legend.text.size = 1.2)
# tmap_save(tm = map_P_26, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/P_map_26.png")
# 
# map_P_85 = tm_shape(World) + tm_fill() + tm_shape(basins_map) + tm_borders() +
#   tm_fill(col="dev_PCP_85", breaks=c(-600,-200,-100,0,100,200,500,10000),
#           labels=c("< -200", "-200 to -100", "-100 to 0", "0 to 100", "100 to 200", "200 to 500", "> 500"),
#           title=expression("Precipitation deviation (mm"~yr^-1*")"),
#           palette="RdBu") +
#   tm_layout(frame=FALSE,legend.position = c("left","bottom"),
#             legend.title.size = 1.3,
#             legend.text.size = 1.2)
# tmap_save(tm = map_P_85, dpi=300,filename = "E:/Papiers/2020_Fabre et al Global carbon futur/P_map_85.png")





bg_col = ifelse(climats$OCEAN == "Arctic Ocean", "blue",
                ifelse(climats$OCEAN == "Atlantic Ocean", "darkgreen",
                       ifelse(climats$OCEAN == "Pacific Ocean", "red",
                              ifelse(climats$OCEAN == "Indian Ocean", "yellow","NULL"))))
xlim_Q = c(-75,220)
xlim_SOC = c(-100,25)
xlim_TMP = c(-15,10)
xlim_PCP = c(-40,60)

ylim_26 = c(-65,300)
ylim_85 = c(-75,800)
# Figure 6
par(mfrow=c(2,4))
# RCP26
plot(climats$diff_Q_26, climats$diff_26, pch=21, xlim=xlim_Q,ylim=ylim_26, col="darkgrey", 
     bg=bg_col, xlab="relative Q change (%)", ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_26~diff_Q_26, data=climats)$coefficients[1] + 
        lm(diff_26~diff_Q_26, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_26~diff_Q_26, data=climats))$r.squared)

plot(climats$diff_SOC_26, climats$diff_26, pch=21, xlim=xlim_SOC,ylim=ylim_26, col="darkgrey", 
     bg=bg_col, xlab="Relative SOC change (%)", ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_26~diff_SOC_26, data=climats)$coefficients[1] + 
        lm(diff_26~diff_SOC_26, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_26~diff_SOC_26, data=climats))$r.squared)

plot(climats$dev_TMP_26, climats$diff_26, pch=21, xlim=xlim_TMP,ylim=ylim_26, col="darkgrey", 
     bg=bg_col, xlab=expression("Temperature deviation ("~Delta*"°C)"), ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_26~dev_TMP_26, data=climats)$coefficients[1] + 
        lm(diff_26~dev_TMP_26, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_26~dev_TMP_26, data=climats))$r.squared)

plot(climats$diff_PCP_26, climats$diff_26, pch=21, xlim=xlim_PCP,ylim=ylim_26, col="darkgrey", 
     bg=bg_col, xlab=expression("Relative precipitation change (%)"), ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_26~diff_PCP_26, data=climats)$coefficients[1] + 
        lm(diff_26~diff_PCP_26, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_26~diff_PCP_26, data=climats))$r.squared)

# RCP85
plot(climats$diff_Q_85, climats$diff_85, pch=21, xlim=xlim_Q,ylim=ylim_85, col="darkgrey", 
     bg=bg_col, xlab="Relative Q change (%)", ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_85~diff_Q_85, data=climats)$coefficients[1] + 
        lm(diff_85~diff_Q_85, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_85~diff_Q_85, data=climats))$r.squared)

plot(climats$diff_SOC_85, climats$diff_85, pch=21, xlim=xlim_SOC,ylim=ylim_85, col="darkgrey", 
     bg=bg_col, xlab="Relative SOC change (%)", ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_85~diff_SOC_85, data=climats)$coefficients[1] + 
        lm(diff_85~diff_SOC_85, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_85~diff_SOC_85, data=climats))$r.squared)

plot(climats$dev_TMP_85, climats$diff_85, pch=21, xlim=xlim_TMP,ylim=ylim_85, col="darkgrey", 
     bg=bg_col, xlab=expression("Temperature deviation ("~Delta*"°C)"), ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_85~dev_TMP_85, data=climats)$coefficients[1] + 
        lm(diff_85~dev_TMP_85, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_85~dev_TMP_85, data=climats))$r.squared)

plot(climats$diff_PCP_85, climats$diff_85, pch=21, xlim=xlim_PCP,ylim=ylim_85, col="darkgrey", 
     bg=bg_col, xlab=expression("Relative precipitation change (%)"), ylab="Relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
lines(x=c(-1000,1000), y=lm(diff_85~diff_PCP_85, data=climats)$coefficients[1] + 
        lm(diff_85~diff_PCP_85, data=climats)$coefficients[2] * c(-1000,1000), lty=2)
print(summary(lm(diff_85~diff_PCP_85, data=climats))$r.squared)







pass_26 = as.data.frame(ddply(climats,~group_lat,summarise,
                              Q_past_gfdl = sum(Q_hist_rcp26_gfdl),
                              Q_past_hadgem = sum(Q_hist_rcp26_hadgem),
                              Q_past_ipsl = sum(Q_hist_rcp26_ipsl),
                              Q_past_miroc = sum(Q_hist_rcp26_miroc),
                              Q_past_noresm = sum(Q_hist_rcp26_noresm)))


DF = data.frame(group_lat=c('[-60,-50)','[-50,-40)','[-40,-30)','[-30,-20)','[-20,-10)','[-10,0)',
                      '[0,10)','[10,20)','[20,30)','[30,40)','[40,50)','[50,60)','[60,70)','[70,80)'),
                latitude=c('50-60S','40-50S','30-40S','20-30S','10-20S','0-10S',
                      '0-10N','10-20N','20-30N','30-40N','40-50N','50-60N','60-70N','70-80N'),
                observations=c(23,47,751,76,264,7256,1474,1297,558,1686,1248,1476,1859,1449))

DF = merge(DF, pass_26, by='group_lat')





# pass_DOC_26 = as.data.frame(ddply(climats,~group,summarise,
#                               DOC_past_gfdl = sum(DOC_hist_rcp26_gfdl_mean),
#                               DOC_past_hadgem = sum(DOC_hist_rcp26_hadgem_mean),
#                               DOC_past_ipsl = sum(DOC_hist_rcp26_ipsl_mean),
#                               DOC_past_miroc = sum(DOC_hist_rcp26_miroc_mean),
#                               DOC_past_noresm = sum(DOC_hist_rcp26_noresm_mean)))

DF_DOC = data.frame(group_lat=c('[-60,-50)','[-50,-40)','[-40,-30)','[-30,-20)','[-20,-10)','[-10,0)',
                        '[0,10)','[10,20)','[20,30)','[30,40)','[40,50)','[50,60)','[60,70)','[70,80)'),
                    latitude=c('50-60S','40-50S','30-40S','20-30S','10-20S','0-10S',
                               '0-10N','10-20N','20-30N','30-40N','40-50N','50-60N','60-70N','70-80N'),
                Fabre_ref=c(0.126565269,0.216620004,3.230961125,0.282813862,1.058627799,39.79689455,
                            8.138175292,7.117522477,2.895710928,7.216698875,5.602196844,7.553932046,
                            10.68751637,4.603530517),
                Fabre_sd=c(0,0.056359776,0.603164072,0.032210837,0.08670067,7.463246529,1.208563863,
                           0.795053409,0.477426282,0.813063442,0.24544418,0.178740929,0.412302479,
                           0.868692268))
df_F_DOC = merge(DF_DOC, DOC_lat, by='group_lat')


# df_F_DOC = as.data.frame(read_xlsx("E:/Papiers/2020_Fabre et al Global carbon futur/New_version_1_model/comp obs mod.xlsx", sheet = 'Sheet1'))
library(reshape2)
df_F_DOC_melt = melt(df_F_DOC)
means = df_F_DOC_melt[df_F_DOC_melt$variable %in% c('Fabre_ref', 'DOC_26_hist'),]
sds = df_F_DOC_melt[df_F_DOC_melt$variable %in% c('Fabre_sd', 'DOC_26_hist_sd'),]
sds[sds$variable=='Fabre_sd',]$variable = 'Fabre_ref'
sds[sds$variable=='DOC_26_hist_sd',]$variable = 'DOC_26_hist'

df_FDOC = merge(means,sds, by=c("group_lat","variable","latitude"))
df_FDOC = df_FDOC[order(df_FDOC$variable), ]

colnames(df_FDOC) = c("group_lat","variable","latitude","value","sd")

df_FDOC$order = c(6,5,4,3,2,1,7,8,9,10,11,12,13,14,6,5,4,3,2,1,7,8,9,10,11,12,13,14)

# Graphique par defaut

p1 = ggplot(df_FDOC, aes(x=reorder(latitude, order), y=value, fill=variable)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(x="Latitude", y = "DOC flux (TgC"~yr^-1*")")+
  theme_publish() +
  scale_fill_manual(name="", values=c('white','darkgrey'),
                    labels=c("Observations","Historical")) +
  theme(legend.position="None")


df_F_Q = as.data.frame(read_xlsx("E:/Papiers/2020_Fabre et al Global carbon futur/New_version_1_model/comp obs mod.xlsx", sheet = 'Feuil1'))



# Graphique par defaut
p2 = ggplot(df_F_Q, aes(x=reorder(latitude, order), y=value, fill=group)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(x="Latitude", y = "Discharge ("~km^3*""~yr^-1*")")+
  theme_publish() +
  scale_fill_manual(name="", values=c('white','darkgrey'),
                    labels=c("Observations","Historical")) +
  theme(legend.justification=c(0,1), legend.position=c(0.7,0.97),
        legend.background = element_rect(fill="white", size=.5))# +
  # geom_text(aes(x=reorder(latitude, order), y=value+sd, label=groups),
  #   position = position_dodge(width = 1), vjust=0.3, hjust = -0.5, size = 3)


ggarrange((p1 | p2) + plot_layout(axis_titles='collect'))







# Test validity of the model ####
data = as.data.frame(read_xlsx("E:/Papiers/2020_Fabre et al Global carbon futur/Test_decades_future_DOC.xlsx", sheet="Feuil2", na="NA"))

final = data %>% group_by(River,Decade) %>%
  summarize(DOC = sum(DOC)) %>% 
  as.data.frame()

row1=c("Garonne","Global",NA,NA,NA,NA,NA,NA,NA)
row2=c("Murray","Global",NA,NA,NA,NA,NA,NA,NA)
row3=c("Solimoes","Global",NA,NA,NA,NA,NA,NA,NA)
row4=c("Yenisei","Global",NA,NA,NA,NA,NA,NA,NA)

final=rbind(final,row1,row2,row3,row4)

final$DOC = NULL
final$K=NA
final$Q=NA
final$Std1=NA
final$Std2=NA
# final$Sigma=NA
final$RMSE=NA
final$R=NA
final$pvalue=NA

list_theta1 = c(3,11,4.5,12)
list_theta2 = c(0.001,0.007,2.5,1)

for (i in 1:length(unique(final$River))){
  df = data[data$River==unique(final$River)[i],]
  print(unique(final$River)[i])
  
  for (j in 1:length(final[final$River==unique(final$River)[i],]$Decade)){
    df_s = df[df$Decade==final[final$River==unique(final$River)[i],]$Decade[j],]
    print(final[final$River==unique(final$River)[i],]$Decade[j])
    print(nrow(df_s))
    river = unique(final$River)[i]
    decade = final[final$River==unique(final$River)[i],]$Decade[j]
    
    start = list(theta1 = list_theta1[i], theta2 = list_theta2[i])
    lower = 0.001
    
    nls1=na.omit(as.matrix(df_s["DOC"]))
    nls2=na.omit(as.matrix(df_s["Q_mm"]))
    formula = (nls1) ~ ((theta1*nls2)/(nls2+theta2))
    
    if(length(nls1)>1){
      nls_DOC=nlxb(formula = formula, control=nls.control(maxiter = 30000, tol = 1e-05),
                   start=start, data=df_s, trace=F, algorithm='port', lower=lower)
      fit.nls <- nls2(formula = formula, data = df_s, start = nls_DOC$coefficients,
                      algorithm = "brute-force")
      
      
      final[final$River==river & final$Decade==decade,]$K=summary(fit.nls)$coeff[1]
      final[final$River==river & final$Decade==decade,]$Q=summary(fit.nls)$coeff[2]
      
      df_s$DOC_pred = ((summary(fit.nls)$coeff[1]*df_s$Q_mm)/(summary(fit.nls)$coeff[2]+df_s$Q_mm))
      
      final[final$River==river & final$Decade==decade,]$R= summary(lm(DOC_pred ~ DOC, df_s))$r.squared
      final[final$River==river & final$Decade==decade,]$RMSE=sqrt(sum(nls1-predict(fit.nls))^2/length(nls1))
      
      final[final$River==river & final$Decade==decade,]$pvalue=(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]+summary(fit.nls)$coeff[,"Pr(>|t|)"][2])/2
      final[final$River==river & final$Decade==decade,]$Std1= summary(fit.nls)$coeff[,"Std. Error"][1]
      if (is.na(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]) == F & summary(fit.nls)$coeff[,"Pr(>|t|)"][1] > 0.95){
        final[final$River==river & final$Decade==decade,]$Std1=0
      }
      final[final$River==river & final$Decade==decade,]$Std2=summary(fit.nls)$coeff[,"Std. Error"][2]
      if (is.na(final[final$River==river & final$Decade==decade,]$Std2) == F & final[final$River==river & final$Decade==decade,]$Std2 > final[final$River==river & final$Decade==decade,]$Q) {
        final[final$River==river & final$Decade==decade,]$Std2 = final[final$River==river & final$Decade==decade,]$Q}
      if (is.na(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]) == F & summary(fit.nls)$coeff[,"Pr(>|t|)"][2] > 0.95){
        final[final$River==river & final$Decade==decade,]$Std2=0
      }
    } else {
      final[final$River==river & final$Decade==decade,]$K=NA
      final[final$River==river & final$Decade==decade,]$Q=NA
      final[final$River==river & final$Decade==decade,]$RMSE=NA
      final[final$River==river & final$Decade==decade,]$R=NA
      final[final$River==river & final$Decade==decade,]$pvalue=NA
      final[final$River==river & final$Decade==decade,]$Std1= NA
      final[final$River==river & final$Decade==decade,]$Std2= NA
    }
  }
}

for (i in 1:length(unique(final$River))){
  df = data[data$River==unique(final$River)[i],]
  print(unique(final$River)[i])
  river = unique(final$River)[i]
  start = list(theta1 = 20, theta2 = 1)
  lower = 0.001
  nls1=na.omit(as.matrix(df["DOC"])) ; nls2=na.omit(as.matrix(df["Q_mm"]))
  formula = (nls1) ~ ((theta1*nls2)/(nls2+theta2))
  
  if(length(nls1)>1){
    nls_DOC=nlxb(formula = formula, control=nls.control(maxiter = 30000, tol = 1e-05),
                 start=start, data=df, trace=F, algorithm='port', lower=lower)
    fit.nls <- nls2(formula = formula, data = df, start = nls_DOC$coefficients,
                    algorithm = "brute-force")
    
    
    final[final$River==river & final$Decade=="Global",]$K=summary(fit.nls)$coeff[1]
    final[final$River==river & final$Decade=="Global",]$Q=summary(fit.nls)$coeff[2]
    
    df$DOC_pred = ((summary(fit.nls)$coeff[1]*df$Q_mm)/(summary(fit.nls)$coeff[2]+df$Q_mm))
    
    final[final$River==river & final$Decade=="Global",]$R= summary(lm(DOC_pred ~ DOC, df))$r.squared
    final[final$River==river & final$Decade=="Global",]$RMSE=sqrt(sum(nls1-predict(fit.nls))^2/length(nls1))
    
    final[final$River==river & final$Decade=="Global",]$pvalue=(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]+summary(fit.nls)$coeff[,"Pr(>|t|)"][2])/2
    final[final$River==river & final$Decade=="Global",]$Std1= summary(fit.nls)$coeff[,"Std. Error"][1]
    if (is.na(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]) == F & summary(fit.nls)$coeff[,"Pr(>|t|)"][1] > 0.95){
      final[final$River==river & final$Decade=="Global",]$Std1=0
    }
    final[final$River==river & final$Decade=="Global",]$Std2=summary(fit.nls)$coeff[,"Std. Error"][2]
    if (is.na(final[final$River==river & final$Decade=="Global",]$Std2) == F & final[final$River==river & final$Decade=="Global",]$Std2 > final[final$River==river & final$Decade=="Global",]$Q) {
      final[final$River==river & final$Decade=="Global",]$Std2 = final[final$River==river & final$Decade=="Global",]$Q}
    if (is.na(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]) == F & summary(fit.nls)$coeff[,"Pr(>|t|)"][2] > 0.95){
      final[final$River==river & final$Decade=="Global",]$Std2=0
    }
  } else {
    final[final$River==river & final$Decade=="Global",]$K=NA
    final[final$River==river & final$Decade=="Global",]$Q=NA
    final[final$River==river & final$Decade=="Global",]$RMSE=NA
    final[final$River==river & final$Decade=="Global",]$R=NA
    final[final$River==river & final$Decade=="Global",]$pvalue=NA
    final[final$River==river & final$Decade=="Global",]$Std1= NA
    final[final$River==river & final$Decade=="Global",]$Std2= NA
  }
}


# Adding column based on other column:
final = final %>%
  mutate(color = case_when(
    Decade == "1970-1979" ~ "darkgrey",
    Decade == "1980-1989" ~ "limegreen",
    Decade == "1990-1999" ~ "blue",
    Decade == "2000-2009" ~ "red",
    Decade == "2010-2019" ~ "orange"),  
    
    shape = case_when(Decade == "1970-1979" ~ 4,
                      Decade == "1980-1989" ~ 15,
                      Decade == "1990-1999" ~ 16,
                      Decade == "2000-2009" ~ 17,
                      Decade == "2010-2019" ~ 18))

final$color[is.na(final$K)] = "transparent"
final$shape[is.na(final$K)] = NA

for (k in 1: length(unique(final$River))){
  
  df = data[data$River==unique(final$River)[k],]
  river = unique(final$River)[k]
  
  x1 = seq(0,1.05*max(na.omit(df$Q_mm)),0.001)
  
  K = final[final$River==river & final$Decade=="1970-1979",]$K
  Q = final[final$River==river & final$Decade=="1970-1979",]$Q
  Std1 = final[final$River==river & final$Decade=="1970-1979",]$Std1
  Std2 = final[final$River==river & final$Decade=="1970-1979",]$Std2
  
  if(is.na(final[final$River==river & final$Decade == "1970-1979",]$K)){
    y0_1970 = rep(0, length(x1)); y1_1970 = rep(0, length(x1)); y2_1970 = rep(0, length(x1))
  } else {
    y0_1970 = (K*seq(0,1.05*max(na.omit(na.omit(df$Q_mm))),0.001))/(Q+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
    y1_1970 = (K+Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q-Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
    y2_1970 = (K-Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q+Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  }
  
  K = final[final$River==river & final$Decade=="1980-1989",]$K
  Q = final[final$River==river & final$Decade=="1980-1989",]$Q
  Std1 = final[final$River==river & final$Decade=="1980-1989",]$Std1
  Std2 = final[final$River==river & final$Decade=="1980-1989",]$Std2
  
  if(is.na(final[final$River==river & final$Decade == "1980-1989",]$K)){
    y0_1980 = rep(0, length(x1)); y1_1980 = rep(0, length(x1)); y2_1980 = rep(0, length(x1))
  } else {
    y0_1980 = (K*seq(0,1.05*max(na.omit(df$Q_mm)),0.001))/(Q+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
    y1_1980 = (K+Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q-Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
    y2_1980 = (K-Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q+Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  }
  
  K = final[final$River==river & final$Decade=="1990-1999",]$K
  Q = final[final$River==river & final$Decade=="1990-1999",]$Q
  Std1 = final[final$River==river & final$Decade=="1990-1999",]$Std1
  Std2 = final[final$River==river & final$Decade=="1990-1999",]$Std2
  
  if(is.na(final[final$River==river & final$Decade == "1990-1999",]$K)){
    y0_1990 = rep(0, length(x1)); y1_1990 = rep(0, length(x1)); y2_1990 = rep(0, length(x1))
  } else {
    y0_1990 = (K*seq(0,1.05*max(na.omit(df$Q_mm)),0.001))/(Q+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
    y1_1990 = (K+Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q-Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
    y2_1990 = (K-Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q+Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  }
  
  K = final[final$River==river & final$Decade=="2000-2009",]$K
  Q = final[final$River==river & final$Decade=="2000-2009",]$Q
  Std1 = final[final$River==river & final$Decade=="2000-2009",]$Std1
  Std2 = final[final$River==river & final$Decade=="2000-2009",]$Std2
  
  y0_2000 = (K*seq(0,1.05*max(na.omit(na.omit(df$Q_mm))),0.001))/(Q+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  y1_2000 = (K+Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q-Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  y2_2000 = (K-Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q+Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  
  K = final[final$River==river & final$Decade=="2010-2019",]$K
  Q = final[final$River==river & final$Decade=="2010-2019",]$Q
  Std1 = final[final$River==river & final$Decade=="2010-2019",]$Std1
  Std2 = final[final$River==river & final$Decade=="2010-2019",]$Std2
  
  y0_2010 = (K*seq(0,1.05*max(na.omit(df$Q_mm)),0.001))/(Q+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  y1_2010 = (K+Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q-Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  y2_2010 = (K-Std1)*seq(0,1.05*max(na.omit(df$Q_mm)),0.001)/(Q+Std2+(seq(0,1.05*max(na.omit(df$Q_mm)),0.001)))
  
  dat = data.frame(x1,
                   y0_1970,y1_1970,y2_1970,
                   y0_1980,y1_1980,y2_1980,
                   y0_1990,y1_1990,y2_1990,
                   y0_2000,y1_2000,y2_2000,
                   y0_2010,y1_2010,y2_2010)
  
  
  p[[k]] = ggplot(df, aes(Q_mm, DOC)) +
    xlab("Daily discharge ("~mm*"."~d^-1*")") +
    ylab("DOC concentration (mg."~L^-1*")") +
    expand_limits(x=0, y=0) +
    ggtitle(paste(river)) +
    geom_line(data=dat, aes(x=x1, y=y0_1970, color=as.factor("1970-1979")), linewidth=1) +
    geom_ribbon(data=dat, aes(x=x1, ymin=y2_1970, ymax=y1_1970), inherit.aes = FALSE, fill="darkgrey", alpha=0.1) +
    geom_line(data=dat, aes(x=x1, y=y0_1980, color=as.factor("1980-1989")), linewidth=1) +
    geom_ribbon(data=dat, aes(x=x1, ymin=y2_1980, ymax=y1_1980), inherit.aes = FALSE, fill="limegreen", alpha=0.1) +
    geom_line(data=dat, aes(x=x1, y=y0_1990, color=as.factor("1990-1999")), linewidth=1) +
    geom_ribbon(data=dat, aes(x=x1, ymin=y2_1990, ymax=y1_1990), inherit.aes = FALSE, fill="blue", alpha=0.1) +
    geom_line(data=dat, aes(x=x1, y=y0_2000, color=as.factor("2000-2009")), linewidth=1) +
    geom_ribbon(data=dat, aes(x=x1, ymin=y2_2000, ymax=y1_2000), inherit.aes = FALSE, fill="red", alpha=0.1) +
    geom_line(data=dat, aes(x=x1, y=y0_2010, color=as.factor("2010-2019")), linewidth=1) +
    geom_ribbon(data=dat, aes(x=x1, ymin=y2_2010, ymax=y1_2010), inherit.aes = FALSE, fill="orange", alpha=0.1) +
    geom_point(data=df, aes(y=DOC, shape=Decade, color=Decade), na.rm = FALSE, size=2.5) +
    theme_pubr() +
    scale_shape_manual(name = "Decade",
                       labels = c("1970-1979","1980-1989","1990-1999","2000-2009","2010-2019"),
                       values = as.vector(final[final$River==unique(final$River)[k],]$shape)) +
    scale_color_manual(name = "Decade",
                       labels = c("1970-1979","1980-1989","1990-1999","2000-2009","2010-2019"),
                       values = as.vector(final[final$River==unique(final$River)[k],]$color)) +
    theme(axis.title.y = element_text(size = 12, angle = 90)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.text = element_text(size = 12))
  
}

# png(paste("E:/Papiers/2020_Fabre et al Global carbon futur/Soumission/Figure_DOC_Q/",river,".png", sep=''), width = 2800, height = 2000, res = 300)

ggarrange(p[[1]], 
          p[[2]], 
          p[[3]], 
          p[[4]], ncol = 2, nrow=2, common.legend = T, legend='bottom')
# dev.off()

write.csv(final, file=paste0(dir,"/final_stats.csv"))





















