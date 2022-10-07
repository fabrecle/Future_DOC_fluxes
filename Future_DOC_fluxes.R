# Introduction ####
library(zoo) ; library(plyr) ; library(dplyr)
library(ggplot2) ; library(ggrepel); library(readxl)
library(stringr) ; library(matrixStats) ; library(lubridate)
library(sf) ; library(raster) ; library(spData) ; library(spDataLarge)

library(tmap) ; library(nlstools)
library(nlmrt) ; library(nls2)


# Figure 3 & Table 1 ####
data = as.data.frame(read_xlsx("Test_decades_future_DOC.xlsx", sheet="Feuil2", na="NA"))

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
final$RMSE=NA
final$R=NA
final$pvalue=NA

for (i in 1:length(unique(final$River))){
  df = data[data$River==unique(final$River)[i],]
  print(unique(final$River)[i])
  
for (j in 1:length(final[final$River==unique(final$River)[i],]$Decade)){
  df_s = df[df$Decade==final[final$River==unique(final$River)[i],]$Decade[j],]
  print(final[final$River==unique(final$River)[i],]$Decade[j])
  
  river = unique(final$River)[i]
  decade = final[final$River==unique(final$River)[i],]$Decade[j]
  
start = list(theta1 = 20, theta2 = 1)
lower = 0.001

nls1=na.omit(as.matrix(df_s["DOC"]))
nls2=na.omit(as.matrix(df_s["Q_mm"]))
formula = (nls1) ~ ((theta1*nls2)/(nls2+theta2))

if(length(nls1)>1){
nls_DOC=nlxb(formula = formula, control=nls.control(maxiter = 1000, tol = 1e-05),
             start=start, data=df_s, trace=F, algorithm='port', lower=lower)
fit.nls <- nls2(formula = formula, data = df_s, start = nls_DOC$coefficients,
                algorithm = "brute-force")


final[final$River==river & final$Decade==decade,]$K=summary(fit.nls)$coeff[1]
final[final$River==river & final$Decade==decade,]$Q=summary(fit.nls)$coeff[2]

df_s$DOC_pred = ((summary(fit.nls)$coeff[1]*df_s$Q_mm)/(summary(fit.nls)$coeff[2]+df_s$Q_mm))

final[final$River==river & final$Decade==decade,]$R= summary(lm(DOC_pred ~ DOC, df_s))$r.squared
final[final$River==river & final$Decade==decade,]$RMSE=sqrt(sum(nls1-predict(fit.nls))^2/length(nls1))
final[final$River==river & final$Decade==decade,]$pvalue=(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]+summary(fit.nls)$coeff[,"Pr(>|t|)"][2])/2
final[final$River==river & final$Decade==decade,]$Std1= summary(fit.nls)$coeff[,"Std. Error"][1]/sqrt(nrow(df_s))
if (is.na(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]) == F & summary(fit.nls)$coeff[,"Pr(>|t|)"][1] > 0.95){
  final[final$River==river & final$Decade==decade,]$Std1=0
}
final[final$River==river & final$Decade==decade,]$Std2=summary(fit.nls)$coeff[,"Std. Error"][2]/sqrt(nrow(df_s))
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
  nls_DOC=nlxb(formula = formula, control=nls.control(maxiter = 1000, tol = 1e-05),
               start=start, data=df, trace=F, algorithm='port', lower=lower)
  fit.nls <- nls2(formula = formula, data = df, start = nls_DOC$coefficients,
                  algorithm = "brute-force")
  
  
  final[final$River==river & final$Decade=="Global",]$K=summary(fit.nls)$coeff[1]
  final[final$River==river & final$Decade=="Global",]$Q=summary(fit.nls)$coeff[2]
  
  df$DOC_pred = ((summary(fit.nls)$coeff[1]*df$Q_mm)/(summary(fit.nls)$coeff[2]+df$Q_mm))
  
  final[final$River==river & final$Decade=="Global",]$R= summary(lm(DOC_pred ~ DOC, df))$r.squared
  final[final$River==river & final$Decade=="Global",]$RMSE=sqrt(sum(nls1-predict(fit.nls))^2/length(nls1))
  
  final[final$River==river & final$Decade=="Global",]$pvalue=(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]+summary(fit.nls)$coeff[,"Pr(>|t|)"][2])/2
  final[final$River==river & final$Decade=="Global",]$Std1= summary(fit.nls)$coeff[,"Std. Error"][1]/sqrt(nrow(df))
  if (is.na(summary(fit.nls)$coeff[,"Pr(>|t|)"][1]) == F & summary(fit.nls)$coeff[,"Pr(>|t|)"][1] > 0.95){
    final[final$River==river & final$Decade=="Global",]$Std1=0
  }
  final[final$River==river & final$Decade=="Global",]$Std2=summary(fit.nls)$coeff[,"Std. Error"][2]/sqrt(nrow(df))
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

# Colors and shapes
for (i in 1:length(final$River)){
if (!is.na(final$K[i])){
  if (final$Decade[i] == "1970-1979"){
    final$color[i] = "darkgrey"; final$shape[i] = 4
  } else if (final$Decade[i] == "1980-1989"){
    final$color[i] = "limegreen"; final$shape[i] = 15
  } else if (final$Decade[i] == "1990-1999"){
    final$color[i] = "blue"; final$shape[i] = 16
  } else if (final$Decade[i] == "2000-2009"){
    final$color[i] = "red"; final$shape[i] = 17
  } else if (final$Decade[i] == "2010-2019"){
    final$color[i] = "orange"; final$shape[i] = 18
  }
} else {
  final$color[i] = "transparent";  final$shape[i] = NA
}
}

for (k in 1: length(unique(final$River))){
  
  df = data[data$River==unique(final$River)[k],]
  river = unique(final$River)[k]

  x1 = seq(0,max(na.omit(df$Q_mm))+0.2,0.001)
  
  K = final[final$River==river & final$Decade=="1970-1979",]$K
  Q = final[final$River==river & final$Decade=="1970-1979",]$Q
  Std1 = final[final$River==river & final$Decade=="1970-1979",]$Std1
  Std2 = final[final$River==river & final$Decade=="1970-1979",]$Std2
  
  if(is.na(final[final$River==river & final$Decade == "1970-1979",]$K)){
    y0_1970 = rep(0, length(x1)); y1_1970 = rep(0, length(x1)); y2_1970 = rep(0, length(x1))
  } else {
    y0_1970 = (K*seq(0,max(na.omit(na.omit(df$Q_mm)))+0.2,0.001))/(Q+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
    y1_1970 = (K+Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q-Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
    y2_1970 = (K-Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q+Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
 }

  K = final[final$River==river & final$Decade=="1980-1989",]$K
  Q = final[final$River==river & final$Decade=="1980-1989",]$Q
  Std1 = final[final$River==river & final$Decade=="1980-1989",]$Std1
  Std2 = final[final$River==river & final$Decade=="1980-1989",]$Std2

  if(is.na(final[final$River==river & final$Decade == "1980-1989",]$K)){
    y0_1980 = rep(0, length(x1)); y1_1980 = rep(0, length(x1)); y2_1980 = rep(0, length(x1))
  } else {
    y0_1980 = (K*seq(0,max(na.omit(df$Q_mm))+0.2,0.001))/(Q+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
    y1_1980 = (K+Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q-Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
    y2_1980 = (K-Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q+Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  }

  K = final[final$River==river & final$Decade=="1990-1999",]$K
  Q = final[final$River==river & final$Decade=="1990-1999",]$Q
  Std1 = final[final$River==river & final$Decade=="1990-1999",]$Std1
  Std2 = final[final$River==river & final$Decade=="1990-1999",]$Std2
    
  if(is.na(final[final$River==river & final$Decade == "1990-1999",]$K)){
    y0_1990 = rep(0, length(x1)); y1_1990 = rep(0, length(x1)); y2_1990 = rep(0, length(x1))
  } else {
    y0_1990 = (K*seq(0,max(na.omit(df$Q_mm))+0.2,0.001))/(Q+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
    y1_1990 = (K+Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q-Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
    y2_1990 = (K-Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q+Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  }

  K = final[final$River==river & final$Decade=="2000-2009",]$K
  Q = final[final$River==river & final$Decade=="2000-2009",]$Q
  Std1 = final[final$River==river & final$Decade=="2000-2009",]$Std1
  Std2 = final[final$River==river & final$Decade=="2000-2009",]$Std2
  
  y0_2000 = (K*seq(0,max(na.omit(na.omit(df$Q_mm)))+0.2,0.001))/(Q+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  y1_2000 = (K+Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q-Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  y2_2000 = (K-Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q+Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  
  K = final[final$River==river & final$Decade=="2010-2019",]$K
  Q = final[final$River==river & final$Decade=="2010-2019",]$Q
  Std1 = final[final$River==river & final$Decade=="2010-2019",]$Std1
  Std2 = final[final$River==river & final$Decade=="2010-2019",]$Std2
  
  y0_2010 = (K*seq(0,max(na.omit(df$Q_mm))+0.2,0.001))/(Q+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  y1_2010 = (K+Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q-Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  y2_2010 = (K-Std1)*seq(0,max(na.omit(df$Q_mm))+0.2,0.001)/(Q+Std2+(seq(0,max(na.omit(df$Q_mm))+0.2,0.001)))
  
  dat = data.frame(x1,
                   y0_1970,y1_1970,y2_1970,
                   y0_1980,y1_1980,y2_1980,
                   y0_1990,y1_1990,y2_1990,
                   y0_2000,y1_2000,y2_2000,
                   y0_2010,y1_2010,y2_2010)


  p = ggplot(df, aes(Q_mm, DOC)) +
    xlab("Daily discharge ("~mm*"."~j^-1*")") +
    ylab("Dissolved organic carbon concentration (mg."~L^-1*")") +
    expand_limits(x=0, y=0) +
    coord_cartesian(xlim = c(0,1.05*max(na.omit(df$Q_mm))), ylim = c(0,1.05*max(na.omit(df$DOC)))) +
    geom_line(data=dat, aes_string(x=x1, y=y0_1970, color=as.factor("1970-1979")), size=1) +
    geom_ribbon(data=dat, aes_string(x=x1, ymin=y2_1970, ymax=y1_1970), inherit.aes = FALSE, fill="darkgrey", alpha=0.3) +
    geom_line(data=dat, aes_string(x=x1, y=y0_1980, color=as.factor("1980-1989")), size=1) +
    geom_ribbon(data=dat, aes_string(x=x1, ymin=y2_1980, ymax=y1_1980), inherit.aes = FALSE, fill="limegreen", alpha=0.3) +
    geom_line(data=dat, aes_string(x=x1, y=y0_1990, color=as.factor("1990-1999")), size=1) +
    geom_ribbon(data=dat, aes_string(x=x1, ymin=y2_1990, ymax=y1_1990), inherit.aes = FALSE, fill="blue", alpha=0.3) +
    geom_line(data=dat, aes_string(x=x1, y=y0_2000, color=as.factor("2000-2009")), size=1) +
    geom_ribbon(data=dat, aes_string(x=x1, ymin=y2_2000, ymax=y1_2000), inherit.aes = FALSE, fill="red", alpha=0.3) +
    geom_line(data=dat, aes_string(x=x1, y=y0_2010, color=as.factor("2010-2019")), size=1) +
    geom_ribbon(data=dat, aes_string(x=x1, ymin=y2_2010, ymax=y1_2010), inherit.aes = FALSE, fill="orange", alpha=0.3) +
    geom_point(data=df, aes(y=DOC, shape=Decade, color=Decade), na.rm = FALSE, size=2.5) +
    theme_minimal() +
    scale_shape_manual(name = "Decade",
                       labels = c("1970-1979","1980-1989","1990-1999","2000-2009","2010-2019"),
                       values = as.vector(final[final$River==unique(final$River)[k],]$shape)) +
    scale_color_manual(name = "Decade",
                       labels = c("1970-1979","1980-1989","1990-1999","2000-2009","2010-2019"),
                       values = as.vector(final[final$River==unique(final$River)[k],]$color)) +
    theme(axis.title.y = element_text(size = 20, angle = 90)) +
    theme(axis.title.x = element_text(size = 20)) +
    theme(axis.text = element_text(size = 15))
}

write.csv(final, file=paste("Table1.csv", sep=""))





# Preparing dataset ####
climats = as.data.frame(read_xlsx("Climats_future_v2.xlsx", sheet="Climats_future"))

rownames(climats) <- climats[,1]
climats[,1] <- NULL

climats = na.omit(climats)

# Importing discharges datasets ####
path = "Files"
path_final = "Final"

path_files = "Future_discharge/Complete_Routing_run_"

list_models = c("gfdl-esm2m_HistAndRcp2p6_1961-2099_v1_25bands",
                "gfdl-esm2m_HistAndRcp8p5_1961-2099_v1_25bands",
                "hadgem2-es_HistAndRcp2p6_1961-2099_25bands",
                "hadgem2-es_HistAndRcp8p5_1961-2099_25bands",
                "ipsl-cm5a-lr_HistAndRcp2p6_1961-2099_v1_25bands",
                "ipsl-cm5a-lr_HistAndRcp8p5_1961-2099_v1_25bands",
                "miroc-esm-chem_HistAndRcp2p6_1961-2099_v1_25bands",
                "miroc-esm-chem_HistAndRcp8p5_1961-2099_v1_25bands",
                "noresm1-m_HistAndRcp2p6_1961-2099_v1_25bands",
                "noresm1-m_HistAndRcp8p5_1961-2099_v1_25bands")

results_Tg_26 = data.frame(river = character(0), climate = character(0),
                           winter_mean_past = double(0), spring_mean_past = double(0),
                           summer_mean_past = double(0), summer_mean_past = double(0),
                           total_mean_past = double(0), total_mean_past = double(0), 
                           total_mean_past = double(0), total_Q_past = double(0),
                           winter_mean_fut = double(0), spring_mean_fut = double(0),
                           summer_mean_fut = double(0), summer_mean_fut = double(0),
                           total_mean_fut = double(0), total_mean_fut = double(0), 
                           total_mean_fut = double(0), total_Q_fut = double(0),
                           av_DOC_past = double(0), av_DOC_fut = double(0))


results_Tg_85 = data.frame(river = character(0), climate = character(0),
                           winter_mean_past = double(0), spring_mean_past = double(0),
                           summer_mean_past = double(0), summer_mean_past = double(0),
                           total_mean_past = double(0), total_mean_past = double(0), 
                           total_mean_past = double(0), total_Q_past = double(0),
                           winter_mean_fut = double(0), spring_mean_fut = double(0),
                           summer_mean_fut = double(0), summer_mean_fut = double(0),
                           total_mean_fut = double(0), total_mean_fut = double(0), 
                           total_mean_fut = double(0), total_Q_fut = double(0),
                           av_DOC_past = double(0), av_DOC_fut = double(0))

# Creating DOC files per watershed : MEAN, MIN and MAX ####
sum = 0; area_tot = 0; iter = 0

for (i in row.names(climats)){
  print(i)
  name = i

  data_watershed = climats[row.names(climats)==i,]
  area = data_watershed$Area

  sum = sum + 1

  #2.6
  table_01 = read.table(file.path(paste(path_files,list_models[1], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_03 = read.table(file.path(paste(path_files,list_models[3], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_05 = read.table(file.path(paste(path_files,list_models[5], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_07 = read.table(file.path(paste(path_files,list_models[7], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_09 = read.table(file.path(paste(path_files,list_models[9], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))

  table_26 = table_01
  table_26$Value = rowMeans(data.frame(table_01$Value,table_05$Value,table_07$Value,table_09$Value))
  names(table_26)[names(table_26) == 'Value'] <- 'Value_26'
  table_26$Std_26 = rowSds(as.matrix(data.frame(table_01$Value,table_05$Value,table_07$Value,table_09$Value)))/5

  rownames(table_26) <- seq(length=nrow(table_26))

  table_26['DOC_conc_26'] = NA
  table_26['DOC_flux_Tg_26_mean'] = NA
  table_26['DOC_flux_Tg_26_min'] = NA
  table_26['DOC_flux_Tg_26_max'] = NA
  
  rows_past = which(table_26$Y <= 2020)
    alpha = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_past))
    beta = (12.6 / ((data_watershed$T_past) + 16.1)) - 0.03
    DOC_conc = alpha * (table_26[rows_past,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / (beta + (table_26[rows_past,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
    table_26[rows_past,'DOC_conc_26'] = DOC_conc
    table_26[rows_past,'DOC_flux_Tg_26_mean'] = DOC_conc *1e-15* table_26[rows_past,'Value_26'] *1000 *60*60*24
    table_26[rows_past,'DOC_flux_Tg_26_min'] = DOC_conc *1e-15* (table_26[rows_past,'Value_26']-table_26[rows_past,'Std_26']) *1000 *60*60*24
    table_26[rows_past,'DOC_flux_Tg_26_max'] = DOC_conc *1e-15* (table_26[rows_past,'Value_26']+table_26[rows_past,'Std_26']) *1000 *60*60*24

  rows_2030 = which(table_26$Y > 2020 & table_26$Y <= 2040)
      alpha = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_past))
      beta = (12.6 / ((data_watershed$T_past) + 16.1)) - 0.03
      DOC_conc = alpha * (table_26[rows_2030,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / (beta + (table_26[rows_2030,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
      table_26[rows_2030,'DOC_conc_26'] = DOC_conc
      table_26[rows_2030,'DOC_flux_Tg_26_mean'] = DOC_conc *1e-15* table_26[rows_2030,'Value_26'] *1000 *60*60*24
      table_26[rows_2030,'DOC_flux_Tg_26_min'] = DOC_conc *1e-15* (table_26[rows_2030,'Value_26']-table_26[rows_2030,'Std_26']) *1000 *60*60*24
      table_26[rows_2030,'DOC_flux_Tg_26_max'] = DOC_conc *1e-15* (table_26[rows_2030,'Value_26']+table_26[rows_2030,'Std_26']) *1000 *60*60*24
      
  rows_2050 = which(table_26$Y > 2040 & table_26$Y <= 2060)
      alpha_mean = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_26_50))
      beta_mean = (12.6 / ((data_watershed$T_26_50) + 16.1)) - 0.03
      DOC_conc = alpha_mean * (table_26[rows_2050,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / ((beta_mean + table_26[rows_2050,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
      table_26[rows_2050,'DOC_conc_26'] = DOC_conc
      table_26[rows_2050,'DOC_flux_Tg_26_mean'] = DOC_conc *1e-15* table_26[rows_2050,'Value_26'] *1000 *60*60*24
      
      alpha_min = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_26_50+data_watershed$P_26_50_std))
      beta_min = (12.6 / ((data_watershed$T_26_50-data_watershed$T_26_50_std) + 16.1)) - 0.03
      DOC_conc = alpha_min * (table_26[rows_2050,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / (beta_min + (table_26[rows_2050,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
      table_26[rows_2050,'DOC_flux_Tg_26_min'] = DOC_conc *1e-15* (table_26[rows_2050,'Value_26']-table_26[rows_2050,'Std_26']) *1000 *60*60*24
      
      alpha_max = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_26_50-data_watershed$P_26_50_std))
      beta_max = (12.6 / ((data_watershed$T_26_50+data_watershed$T_26_50_std) + 16.1)) - 0.03
      DOC_conc = alpha_max * (table_26[rows_2050,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / (beta_max + (table_26[rows_2050,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
      table_26[rows_2050,'DOC_flux_Tg_26_max'] = DOC_conc *1e-15* (table_26[rows_2050,'Value_26']+table_26[rows_2050,'Std_26']) *1000 *60*60*24
      
  rows_2070 = which(table_26$Y > 2060)
    alpha_mean = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_26_70))
    beta_mean = (12.6 / ((data_watershed$T_26_70) + 16.1)) - 0.03
    DOC_conc = alpha_mean * (table_26[rows_2070,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / (beta_mean + (table_26[rows_2070,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
    table_26[rows_2070,'DOC_conc_26'] = DOC_conc
    table_26[rows_2070,'DOC_flux_Tg_26_mean'] = DOC_conc *1e-15* table_26[rows_2070,'Value_26'] *1000 *60*60*24
    
    alpha_min = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_26_70+data_watershed$P_26_70_std))
    beta_min = (12.6 / ((data_watershed$T_26_70-data_watershed$T_26_70_std) + 16.1)) - 0.03
    DOC_conc = alpha_min * (table_26[rows_2070,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / (beta_min + (table_26[rows_2070,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
    table_26[rows_2070,'DOC_flux_Tg_26_min'] = DOC_conc *1e-15* (table_26[rows_2070,'Value_26']-table_26[rows_2070,'Std_26']) *1000 *60*60*24
    
    alpha_max = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_26_70-data_watershed$P_26_70_std))
    beta_max = (12.6 / ((data_watershed$T_26_70+data_watershed$T_26_70_std) + 16.1)) - 0.03
    DOC_conc = alpha_max * (table_26[rows_2070,'Value_26']*86400*1000/(data_watershed$Area*1e6)) / (beta_max + (table_26[rows_2070,'Value_26']*86400*1000/(data_watershed$Area*1e6)))
    table_26[rows_2070,'DOC_flux_Tg_26_max'] = DOC_conc *1e-15* (table_26[rows_2070,'Value_26']+table_26[rows_2070,'Std_26']) *1000 *60*60*24

  #8.5
  table_02 = read.table(file.path(paste(path_files,list_models[2], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_04 = read.table(file.path(paste(path_files,list_models[4], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_06 = read.table(file.path(paste(path_files,list_models[6], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_08 = read.table(file.path(paste(path_files,list_models[8], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))
  table_10 = read.table(file.path(paste(path_files,list_models[10], sep=""),data_watershed$future_file), header=F, sep="", col.names = c("Y","M","D","Value"))

  table_85 = table_02
  table_85$Value = rowMeans(data.frame(table_02$Value,table_04$Value,table_06$Value,table_08$Value,table_10$Value))
  names(table_85)[names(table_85) == 'Value'] <- 'Value_85'
  table_85$Std_85 = rowSds(as.matrix(data.frame(table_02$Value,table_04$Value,table_06$Value,table_08$Value,table_10$Value)))/5
  
  rownames(table_85) <- seq(length=nrow(table_85))
  
  table_85['DOC_conc_85'] = NA
  table_85['DOC_flux_Tg_85_mean'] = NA
  table_85['DOC_flux_Tg_85_min'] = NA
  table_85['DOC_flux_Tg_85_max'] = NA
  
  rows_past = which(table_85$Y <= 2020)
  alpha = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_past))
  beta = (12.6 / ((data_watershed$T_past) + 16.1)) - 0.03
  
  DOC_conc = alpha * (table_85[rows_past,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta + (table_85[rows_past,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_past,'DOC_conc_85'] = DOC_conc
  table_85[rows_past,'DOC_flux_Tg_85_mean'] = DOC_conc *1e-15* table_85[rows_past,'Value_85'] *1000 *60*60*24
  table_85[rows_past,'DOC_flux_Tg_85_min'] = DOC_conc *1e-15* (table_85[rows_past,'Value_85']-table_85[rows_past,'Std_85']) *1000 *60*60*24
  table_85[rows_past,'DOC_flux_Tg_85_max'] = DOC_conc *1e-15* (table_85[rows_past,'Value_85']+table_85[rows_past,'Std_85']) *1000 *60*60*24
  
  rows_2030 = which(table_85$Y > 2020 & table_85$Y <= 2040)
  alpha = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_past))
  beta = (12.6 / ((data_watershed$T_past) + 16.1)) - 0.03
  DOC_conc = alpha * (table_85[rows_2030,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta + (table_85[rows_2030,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_2030,'DOC_conc_85'] = DOC_conc
  table_85[rows_2030,'DOC_flux_Tg_85_mean'] = DOC_conc *1e-15* table_85[rows_2030,'Value_85'] *1000 *60*60*24
  table_85[rows_2030,'DOC_flux_Tg_85_min'] = DOC_conc *1e-15* (table_85[rows_2030,'Value_85']-table_85[rows_2030,'Std_85']) *1000 *60*60*24
  table_85[rows_2030,'DOC_flux_Tg_85_max'] = DOC_conc *1e-15* (table_85[rows_2030,'Value_85']+table_85[rows_2030,'Std_85']) *1000 *60*60*24
  
  rows_2050 = which(table_85$Y > 2040 & table_85$Y <= 2060)
  alpha_mean = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_85_50))
  beta_mean = (12.6 / ((data_watershed$T_85_50) + 16.1)) - 0.03
  DOC_conc = alpha_mean * (table_85[rows_2050,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta_mean + (table_85[rows_2050,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_2050,'DOC_conc_85'] = DOC_conc
  table_85[rows_2050,'DOC_flux_Tg_85_mean'] = DOC_conc *1e-15* table_85[rows_2050,'Value_85'] *1000 *60*60*24
  
  alpha_min = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_85_50+data_watershed$P_85_50_std))
  beta_min = (12.6 / ((data_watershed$T_85_50-data_watershed$T_85_50_std) + 16.1)) - 0.03
  DOC_conc = alpha_min * (table_85[rows_2050,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta_min + (table_85[rows_2050,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_2050,'DOC_flux_Tg_85_min'] = DOC_conc *1e-15* (table_85[rows_2050,'Value_85']-table_85[rows_2050,'Std_85']) *1000 *60*60*24

  alpha_max = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_85_50-data_watershed$P_85_50_std))
  beta_max = (12.6 / ((data_watershed$T_85_50+data_watershed$T_85_50_std) + 16.1)) - 0.03
  DOC_conc = alpha_max * (table_85[rows_2050,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta_max + (table_85[rows_2050,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_2050,'DOC_flux_Tg_85_max'] = DOC_conc *1e-15* (table_85[rows_2050,'Value_85']+table_85[rows_2050,'Std_85']) *1000 *60*60*24
  
  rows_2070 = which(table_85$Y > 2060)
  alpha_mean = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_85_70))
  beta_mean = (12.6 / ((data_watershed$T_85_70) + 16.1)) - 0.03
  DOC_conc = alpha_mean * (table_85[rows_2070,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta_mean + (table_85[rows_2070,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_2070,'DOC_conc_85'] = DOC_conc
  table_85[rows_2070,'DOC_flux_Tg_85_mean'] = DOC_conc *1e-15* table_85[rows_2070,'Value_85'] *1000 *60*60*24
  
  alpha_min = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_85_70+data_watershed$P_85_70_std))
  beta_min = (12.6 / ((data_watershed$T_85_70-data_watershed$T_85_70_std) + 16.1)) - 0.03
  DOC_conc = alpha_min * (table_85[rows_2070,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta_min + (table_85[rows_2070,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_2070,'DOC_flux_Tg_85_min'] = DOC_conc *1e-15* (table_85[rows_2070,'Value_85']-table_85[rows_2070,'Std_85']) *1000 *60*60*24
  
  alpha_max = 4.34+5.21*((data_watershed$soc_past)/(data_watershed$P_85_70-data_watershed$P_85_70_std))
  beta_max = (12.6 / ((data_watershed$T_85_70+data_watershed$T_85_70_std) + 16.1)) - 0.03
  DOC_conc = alpha_max * (table_85[rows_2070,'Value_85']*86400*1000/(data_watershed$Area*1e6)) / (beta_max + (table_85[rows_2070,'Value_85']*86400*1000/(data_watershed$Area*1e6)))
  table_85[rows_2070,'DOC_flux_Tg_85_max'] = DOC_conc *1e-15* (table_85[rows_2070,'Value_85']+table_85[rows_2070,'Std_85']) *1000 *60*60*24

  table_merge <- as.data.frame(merge(table_26, table_85, by= c('Y','M','D')))
 
  write.csv(table_merge, file=file.path(path_final,paste(row.names(data_watershed),".csv", sep='')))
  
  
  
  # summarise scenario 2.6 for past and future periods
  table_26$date <- as.Date(with(table_26, paste(Y, M, D,sep="-")), "%Y-%m-%d")
  yq <- as.yearqtr(as.yearmon(table_26$date, "%Y-%m-%d"))
  table_26$season <- factor(format(yq, "%q"), levels = 1:4,
                      labels = c("winter", "spring", "summer", "fall"))
  table_26$julian <- as.POSIXlt(table_26$date)$yday + 1

  table_26_past <- table_26[year(table_26$date) >= 1971 & year(table_26$date) <= 2000,]
  table_26_fut <- table_26[year(table_26$date) >= 2071,]
  
  table_group_pass_past = as.data.frame(ddply(table_26_past,~julian~season,summarise,
                                              av_DOC_pass = mean(DOC_conc_26, na.rm=T),
                                              DOC_flux_Tg_26_mean=mean(DOC_flux_Tg_26_mean, na.rm=T),
                                              DOC_flux_Tg_26_min =mean(DOC_flux_Tg_26_min, na.rm=T),
                                              DOC_flux_Tg_26_max =mean(DOC_flux_Tg_26_max, na.rm=T),
                                              Value = mean(Value_26, na.rm=T)*86400/1e9,
                                              Value_min = mean(Value_26-Std_26, na.rm=T)*86400/1e9,
                                              Value_max = mean(Value_26+Std_26, na.rm=T)*86400/1e9))
  
  table_group_past = as.data.frame(ddply(table_group_pass_past,~season,summarise,
                                         DOC_flux_Tg_26_mean=sum(DOC_flux_Tg_26_mean, na.rm=T),
                                         DOC_flux_Tg_26_min =sum(DOC_flux_Tg_26_min, na.rm=T),
                                         DOC_flux_Tg_26_max =sum(DOC_flux_Tg_26_max, na.rm=T),
                                         Q = sum(Value, na.rm=T),
                                         Q_min = sum(Value_min, na.rm=T),
                                         Q_max = sum(Value_max, na.rm=T),
                                         av_DOC = mean(av_DOC_pass, na.rm=T)))
  
  table_group_pass_fut = as.data.frame(ddply(table_26_fut,~julian~season,summarise,
                                             av_DOC_pass = mean(DOC_conc_26, na.rm=T),
                                             DOC_flux_Tg_26_mean=mean(DOC_flux_Tg_26_mean, na.rm=T),
                                             DOC_flux_Tg_26_min =mean(DOC_flux_Tg_26_min, na.rm=T),
                                             DOC_flux_Tg_26_max =mean(DOC_flux_Tg_26_max, na.rm=T),
                                             Value = mean(Value_26, na.rm=T)*86400/1e9,
                                             Value_min = mean(Value_26-Std_26, na.rm=T)*86400/1e9,
                                             Value_max = mean(Value_26+Std_26, na.rm=T)*86400/1e9))
  table_group_fut = as.data.frame(ddply(table_group_pass_fut,~season,summarise,
                                        DOC_flux_Tg_26_mean=sum(DOC_flux_Tg_26_mean, na.rm=T),
                                         DOC_flux_Tg_26_min =sum(DOC_flux_Tg_26_min, na.rm=T),
                                         DOC_flux_Tg_26_max =sum(DOC_flux_Tg_26_max, na.rm=T),
                                        Q = sum(Value, na.rm=T),
                                        Q_min = sum(Value_min, na.rm=T),
                                        Q_max = sum(Value_max, na.rm=T),
                                        av_DOC = mean(av_DOC_pass, na.rm=T)))
  
  myList_Tg_26 =list(river=name, climate = climats[row.names(climats)==name,][['CLIMAT']],
                   winter_mean_past = table_group_past[table_group_past[,1]=='winter',][[2]],
                   spring_mean_past = table_group_past[table_group_past[,1]=='spring',][[2]],
                   summer_mean_past = table_group_past[table_group_past[,1]=='summer',][[2]],
                   fall_mean_past = table_group_past[table_group_past[,1]=='fall',][[2]],
                   total_mean_past = sum(table_group_past[[2]]),
                   total_min_past = sum(table_group_past[[3]]),
                   total_max_past = sum(table_group_past[[4]]),
                   total_Q_past = sum(table_group_past[[5]]),
                   total_Q_min_past = sum(table_group_past[[6]]),
                   total_Q_max_past = sum(table_group_past[[7]]),
                   av_DOC_past = mean(table_group_past[[8]]),
                   winter_mean_fut = table_group_fut[table_group_fut[,1]=='winter',][[2]],
                   spring_mean_fut = table_group_fut[table_group_fut[,1]=='spring',][[2]],
                   summer_mean_fut = table_group_fut[table_group_fut[,1]=='summer',][[2]],
                   fall_mean_fut = table_group_fut[table_group_fut[,1]=='fall',][[2]],
                   total_mean_fut = sum(table_group_fut[[2]]),
                   total_min_fut = sum(table_group_fut[[3]]),
                   total_max_fut = sum(table_group_fut[[4]]),
                   total_Q_fut = sum(table_group_fut[[5]]),
                   total_Q_min_fut = sum(table_group_fut[[6]]),
                   total_Q_max_fut = sum(table_group_fut[[7]]),
                   av_DOC_fut = mean(table_group_fut[[8]]))
  results_Tg_26 = rbind(results_Tg_26,myList_Tg_26)
  results_Tg_26$river = as.character(results_Tg_26$river)
  results_Tg_26$climate = as.character(results_Tg_26$climate)

  # summarise scenario 8.5
  table_85$date <- as.Date(with(table_85, paste(Y, M, D,sep="-")), "%Y-%m-%d")
  yq <- as.yearqtr(as.yearmon(table_85$date, "%Y-%m-%d"))
  table_85$season <- factor(format(yq, "%q"), levels = 1:4,
                             labels = c("winter", "spring", "summer", "fall"))
  table_85$julian <- as.POSIXlt(table_85$date)$yday + 1

  table_85_past <- table_85[year(table_85$date) >= 1971 & year(table_85$date) <= 2000,]
  table_85_fut <- table_85[year(table_85$date) >= 2071,]
  
  table_group_pass_past = as.data.frame(ddply(table_85_past,~julian~season,summarise,
                                              av_DOC_pass = mean(DOC_conc_85, na.rm=T),
                                              DOC_flux_Tg_85_mean=mean(DOC_flux_Tg_85_mean, na.rm=T),
                                              DOC_flux_Tg_85_min =mean(DOC_flux_Tg_85_min, na.rm=T),
                                              DOC_flux_Tg_85_max =mean(DOC_flux_Tg_85_max, na.rm=T),
                                              Value = mean(Value_85, na.rm=T)*86400/1e9,
                                              Value_min = mean(Value_85-Std_85, na.rm=T)*86400/1e9,
                                              Value_max = mean(Value_85+Std_85, na.rm=T)*86400/1e9))
  table_group_past = as.data.frame(ddply(table_group_pass_past,~season,summarise,
                                         DOC_flux_Tg_85_mean=sum(DOC_flux_Tg_85_mean, na.rm=T),
                                         DOC_flux_Tg_85_min =sum(DOC_flux_Tg_85_min, na.rm=T),
                                         DOC_flux_Tg_85_max =sum(DOC_flux_Tg_85_max, na.rm=T),
                                         Q = sum(Value, na.rm=T),
                                         Q_min = sum(Value_min, na.rm=T),
                                         Q_max = sum(Value_max, na.rm=T),
                                         av_DOC = mean(av_DOC_pass, na.rm=T)))
  
  table_group_pass_fut = as.data.frame(ddply(table_85_fut,~julian~season,summarise,
                                             av_DOC_pass = mean(DOC_conc_85, na.rm=T),
                                             DOC_flux_Tg_85_mean=mean(DOC_flux_Tg_85_mean, na.rm=T),
                                             DOC_flux_Tg_85_min =mean(DOC_flux_Tg_85_min, na.rm=T),
                                             DOC_flux_Tg_85_max =mean(DOC_flux_Tg_85_max, na.rm=T),
                                             Value = mean(Value_85, na.rm=T)*86400/1e9,
                                             Value_min = mean(Value_85-Std_85, na.rm=T)*86400/1e9,
                                             Value_max = mean(Value_85+Std_85, na.rm=T)*86400/1e9))
  table_group_fut = as.data.frame(ddply(table_group_pass_fut,~season,summarise,DOC_flux_Tg_85_mean=sum(DOC_flux_Tg_85_mean, na.rm=T),
                                        DOC_flux_Tg_85_min =sum(DOC_flux_Tg_85_min, na.rm=T),
                                        DOC_flux_Tg_85_max =sum(DOC_flux_Tg_85_max, na.rm=T),
                                        Q = sum(Value, na.rm=T),
                                        Q_min = sum(Value_min, na.rm=T),
                                        Q_max = sum(Value_max, na.rm=T),
                                        av_DOC = mean(av_DOC_pass, na.rm=T)))
  
  myList_Tg_85=list(river=name, climate = climats[row.names(climats)==name,][['CLIMAT']],
                    latitude = climats[row.names(climats)==name,][['LAT_mouth']],
                    winter_mean_past = table_group_past[table_group_past[,1]=='winter',][[2]],
                    spring_mean_past = table_group_past[table_group_past[,1]=='spring',][[2]],
                    summer_mean_past = table_group_past[table_group_past[,1]=='summer',][[2]],
                    fall_mean_past = table_group_past[table_group_past[,1]=='fall',][[2]],
                    total_mean_past = sum(table_group_past[[2]]),
                    total_min_past = sum(table_group_past[[3]]),
                    total_max_past = sum(table_group_past[[4]]),
                    total_Q_past = sum(table_group_past[[5]]),
                    total_Q_min_past = sum(table_group_past[[6]]),
                    total_Q_max_past = sum(table_group_past[[7]]),
                    av_DOC_past = mean(table_group_past[[8]]),
                    winter_mean_fut = table_group_fut[table_group_fut[,1]=='winter',][[2]],
                    spring_mean_fut = table_group_fut[table_group_fut[,1]=='spring',][[2]],
                    summer_mean_fut = table_group_fut[table_group_fut[,1]=='summer',][[2]],
                    fall_mean_fut = table_group_fut[table_group_fut[,1]=='fall',][[2]],
                    total_mean_fut = sum(table_group_fut[[2]]),
                    total_min_fut = sum(table_group_fut[[3]]),
                    total_max_fut = sum(table_group_fut[[4]]),
                    total_Q_fut = sum(table_group_fut[[5]]),
                    total_Q_min_fut = sum(table_group_fut[[6]]),
                    total_Q_max_fut = sum(table_group_fut[[7]]),
                    av_DOC_fut = mean(table_group_fut[[8]]))
  results_Tg_85 = rbind(results_Tg_85,myList_Tg_85)
  results_Tg_85$river = as.character(results_Tg_85$river)
  results_Tg_85$climate = as.character(results_Tg_85$climate)
}

write.csv(results_Tg_26, file=paste("flux_DOC_by_rivers_Tg_26.csv", sep=""))
write.csv(results_Tg_85, file=paste("flux_DOC_by_rivers_Tg_85.csv", sep=""))


# Figure 7####
plot(100*(results_Tg_26$total_Q_fut - results_Tg_26$total_Q_past)/results_Tg_26$total_Q_past,
     100*(results_Tg_26$total_mean_fut - results_Tg_26$total_mean_past)/results_Tg_26$total_mean_past,
     pch=21,
     xlim=c(-100,150),ylim=c(-100,150), col="darkgrey", bg=ifelse(results_Tg_85$climate == "A", "darkgreen", 
                                                                  ifelse(results_Tg_85$climate == "B", "yellow", 
                                                                         ifelse(results_Tg_85$climate == "C", "limegreen",
                                                                                ifelse(results_Tg_85$climate == "D", "blue",
                                                                                       ifelse(results_Tg_85$climate == "E", "cyan","NULL"))))),
     xlab="relative Q change (%)", ylab="relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
legend(x = "bottomright",          # Position
       legend = c("Tropical", "Semi-arid", "Temperate", "Cold", "Polar"),  # Legend texts
       pch = c(21,21,21,21,21),           # Line types
       col="darkgrey",
       pt.bg = c("darkgreen", "yellow", "limegreen", "blue", "cyan"))

plot(100*(results_Tg_26$av_DOC_fut - results_Tg_26$av_DOC_past)/results_Tg_26$av_DOC_past,
     100*(results_Tg_26$total_mean_fut - results_Tg_26$total_mean_past)/results_Tg_26$total_mean_past,
     pch=21,
     xlim=c(-100,150),ylim=c(-100,150), col="darkgrey", bg=ifelse(results_Tg_85$climate == "A", "darkgreen", 
                                                                  ifelse(results_Tg_85$climate == "B", "yellow", 
                                                                         ifelse(results_Tg_85$climate == "C", "limegreen",
                                                                                ifelse(results_Tg_85$climate == "D", "blue",
                                                                                       ifelse(results_Tg_85$climate == "E", "cyan","NULL"))))),
       xlab="relative DOC conc. change (%)", ylab="relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
legend(x = "bottomright",          # Position
       legend = c("Tropical", "Semi-arid", "Temperate", "Cold", "Polar"),  # Legend texts
       pch = c(21,21,21,21,21),           # Line types
       col="darkgrey",
       pt.bg = c("darkgreen", "yellow", "limegreen", "blue", "cyan"))



plot(100*(results_Tg_85$total_Q_fut - results_Tg_85$total_Q_past)/results_Tg_85$total_Q_past,
     100*(results_Tg_85$total_mean_fut - results_Tg_85$total_mean_past)/results_Tg_85$total_mean_past,
     pch=21,
     xlim=c(-100,350),ylim=c(-100,700), col="darkgrey", bg=ifelse(results_Tg_85$climate == "A", "darkgreen", 
                                                                 ifelse(results_Tg_85$climate == "B", "yellow", 
                                                                        ifelse(results_Tg_85$climate == "C", "limegreen",
                                                                               ifelse(results_Tg_85$climate == "D", "blue",
                                                                                      ifelse(results_Tg_85$climate == "E", "cyan","NULL"))))),
     xlab="relative Q change (%)", ylab="relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
legend(x = "bottomright",          # Position
       legend = c("Tropical", "Semi-arid", "Temperate", "Cold", "Polar"),  # Legend texts
       pch = c(21,21,21,21,21),           # Line types
       col="darkgrey",
       pt.bg = c("darkgreen", "yellow", "limegreen", "blue", "cyan"))




plot(100*(results_Tg_85$av_DOC_fut - results_Tg_85$av_DOC_past)/results_Tg_85$av_DOC_past,
     100*(results_Tg_85$total_mean_fut - results_Tg_85$total_mean_past)/results_Tg_85$total_mean_past,
     pch=21,
     xlim=c(-100,350),ylim=c(-100,700), col="darkgrey", bg=ifelse(results_Tg_85$climate == "A", "darkgreen", 
                                                   ifelse(results_Tg_85$climate == "B", "yellow", 
                                                          ifelse(results_Tg_85$climate == "C", "limegreen",
                                                                 ifelse(results_Tg_85$climate == "D", "blue",
                                                                        ifelse(results_Tg_85$climate == "E", "cyan","NULL"))))),
     xlab="relative DOC conc. change (%)", ylab="relative DOC flux change (%)")
lines(x=c(0,0),y=c(-1000,1000))
lines(x=c(-1000,1000),y=c(0,0))
legend(x = "bottomright",          # Position
       legend = c("Tropical", "Semi-arid", "Temperate", "Cold", "Polar"),  # Legend texts
       pch = c(21,21,21,21,21),           # Line types
       col="darkgrey",
       pt.bg = c("darkgreen", "yellow", "limegreen", "blue", "cyan"))



# Figure 4 ####
results_Tg_26$group = cut(results_Tg_26$latitude, breaks=c(-60, -50, -40, -30, -20, -10, 0,
                                                           10, 20, 30, 40, 50, 60, 70, 80), right = FALSE)
results_Tg_85$group = cut(results_Tg_85$latitude, breaks=c(-60, -50, -40, -30, -20, -10, 0,
                                                           10, 20, 30, 40, 50, 60, 70, 80), right = FALSE)

write.csv(results_Tg_26, file=paste("flux_DOC_by_rivers_Tg_26.csv", sep=""))
write.csv(results_Tg_85, file=paste("flux_DOC_by_rivers_Tg_85.csv", sep=""))

pass_26 = as.data.frame(ddply(results_Tg_26,~group,summarise,
                              Q_past_26 = sum(total_Q_past),
                              Q_min_past_26 = sum(total_Q_min_past),
                              Q_max_past_26 = sum(total_Q_max_past)))

pass_85 = as.data.frame(ddply(results_Tg_85,~group,summarise,
                    Q_past_85 = sum(total_Q_past),
                    Q_min_past_85 = sum(total_Q_min_past),
                    Q_max_past_85 = sum(total_Q_max_past)))


DF = data.frame(group=c('[-60,-50)','[-50,-40)','[-40,-30)','[-30,-20)','[-20,-10)','[-10,0)',
                      '[0,10)','[10,20)','[20,30)','[30,40)','[40,50)','[50,60)','[60,70)','[70,80)'),
                latitude=c('50-60S','40-50S','30-40S','20-30S','10-20S','0-10S',
                      '0-10N','10-20N','20-30N','30-40N','40-50N','50-60N','60-70N','70-80N'),
                observations=c(23,47,751,76,264,7256,1474,1297,558,1686,1248,1476,1859,1449))


DF = merge(DF, pass_26, by='group')
DF = merge(DF, pass_85, by='group')


pass_DOC_26 = as.data.frame(ddply(results_Tg_26,~group,summarise,
                                  DOC_past_26 = sum(total_mean_past),
                                  DOC_min_past_26 = sum(total_min_past),
                                  DOC_max_past_26 = sum(total_max_past)))
pass_DOC_85 = as.data.frame(ddply(results_Tg_85,~group,summarise,
                                  DOC_past_85 = sum(total_mean_past),
                                  DOC_min_past_85 = sum(total_min_past),
                                  DOC_max_past_85 = sum(total_max_past)))
  
df_F_DOC = as.data.frame(read_xlsx("DOC_comp_v2.xlsx"))

df_F_DOC$groups <- c("a","a","a","a","a","a","a","a","a","a","a","a","a","a",
                "b","a","a","a","b","a","a","a","a","a","b","b","b","b",
                "b","a","a","a","b","a","a","a","a","a","b","b","b","b")


ggplot(df_F_DOC, aes(x=reorder(latitude, order), y=value, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(x="Latitude", y = "DOC flux (TgC"~yr^-1*")")+
  theme_classic() +
  scale_fill_manual(name="", values=c('white','lightgrey','darkgrey'),
                    labels=c("Fabre et al. (2020)","Control RCP 2.6","Control RCP 8.5")) +
  theme(legend.justification=c(0,1), legend.position=c(0.7,0.97),
        legend.background = element_rect(fill="white", size=.5)) +
  geom_text(aes(x=reorder(latitude, order), y=value+sd, label=groups),
            position = position_dodge(width = 1), vjust=0.3, hjust = -0.5, size = 3)

df_F_DOC_2 = df_F_DOC %>% group_by(group)  %>%
  summarise(value = sum(value),
            sd = sqrt(sum(sd^2)))

df_F_DOC_2$groups <- c("a","a","a")

ggplot(df_F_DOC_2, aes(x=group, y=value, fill=group)) + 
  ggtitle("Total DOC export (TgC"~yr^-1*")") +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), show.legend = FALSE) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  ylim(0,110) +
  coord_flip() +
  labs(x = "", y = "")+
  theme_classic() +
  scale_fill_manual(name="", values=c('white','lightgrey','darkgrey')) +
  theme(legend.justification=c(0,1), legend.position=c(0.75,0.3),
        legend.background = element_rect(fill="white", size=.5),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=group, y=value+sd, label=groups),
            position = position_dodge(width = 1), vjust=0.3, hjust = -0.5, size = 3)


df2 = as.data.frame(read_xlsx("comparison_obs_model.xlsx"))

df2$groups <- c("a","a","a","a","a","a","a","a","a","a","a","a","a","a",
                "b","b","b","b","b","a","b","b","b","b","b","b","a","b",
                "b","b","b","b","b","a","b","b","b","b","b","b","a","b")

ggplot(df2, aes(x=reorder(latitude, order), y=value, fill=group)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  coord_flip() +
  labs(x="Latitude", y = "Discharge ("~km^3*""~yr^-1*")")+
  theme_classic() +
  scale_fill_manual(name="", values=c('white','lightgrey','darkgrey'),
                    labels=c("Observations","Control RCP 2.6","Control RCP 8.5")) +
  theme(legend.justification=c(0,1), legend.position=c(0.7,0.97),
        legend.background = element_rect(fill="white", size=.5)) +
  geom_text(aes(x=reorder(latitude, order), y=value+sd, label=groups),
    position = position_dodge(width = 1), vjust=0.3, hjust = -0.5, size = 3)

#Full in one bar
df3 = df2 %>% group_by(group)  %>%
  summarise(value = sum(value),
            sd = sqrt(sum(sd^2)))

df3$groups <- c("a","b","b")

ggplot(df3, aes(x=group, y=value, fill=group)) + 
  ggtitle("Total discharge ("~km^3*""~yr^-1*")") +
  geom_bar(stat="identity", color="black", 
           position=position_dodge(), show.legend = FALSE) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) +
  ylim(0,26000) +
  coord_flip() +
  labs(x = "", y = "")+
  theme_classic() +
  scale_fill_manual(name="", values=c('white','lightgrey','darkgrey')) +
  theme(legend.justification=c(0,1), legend.position=c(0.75,0.3),
        legend.background = element_rect(fill="white", size=.5),
        axis.text.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(x=group, y=value+sd, label=groups),
            position = position_dodge(width = 1), vjust=0.3, hjust = -0.5, size = 3)





# Calculate global and regional fluxes ####
files_climat <- as.character(list.files(path=path_final))
names = as.character(strsplit(files_climat,".csv"))

# Global flux
global_Tg = data.frame(date= seq(as.Date("1961/1/1"), as.Date("2099/12/31"), "days"),
                       global_Tg_26_mean = rep(0, 50769),
                       global_Tg_26_min = rep(0, 50769),
                       global_Tg_26_max = rep(0, 50769),
                       global_Tg_85_mean = rep(0, 50769),
                       global_Tg_85_min = rep(0, 50769),
                       global_Tg_85_max = rep(0, 50769))

# Distribution by climate
subset_A = subset(row.names(climats), climats$CLIMAT == 'A' & climats$OCEAN != 'Caspian Sea' & is.na(climats$OCEAN) == F) 
subset_B = subset(row.names(climats), climats$CLIMAT == 'B' & climats$OCEAN != 'Caspian Sea' & is.na(climats$OCEAN) == F) 
subset_C = subset(row.names(climats), climats$CLIMAT == 'C' & climats$OCEAN != 'Caspian Sea' & is.na(climats$OCEAN) == F) 
subset_D = subset(row.names(climats), climats$CLIMAT == 'D' & climats$OCEAN != 'Caspian Sea' & is.na(climats$OCEAN) == F) 
subset_E = subset(row.names(climats), climats$CLIMAT == 'E' & climats$OCEAN != 'Caspian Sea' & is.na(climats$OCEAN) == F) 

area_A = 0 ; area_B = 0 ; area_C = 0 ; area_D = 0 ; area_E = 0
n_A = 0; n_B = 0; n_C = 0; n_D = 0; n_E = 0

climats_Tg = data.frame(date= seq(as.Date("1961/1/1"), as.Date("2099/12/31"), "days"),
                        climat_A_26_mean = rep(0, 50769),
                        climat_A_26_min = rep(0, 50769),
                        climat_A_26_max = rep(0, 50769),
                        climat_A_85_mean = rep(0, 50769),
                        climat_A_85_min = rep(0, 50769),
                        climat_A_85_max = rep(0, 50769),
                        climat_B_26_mean = rep(0, 50769),
                        climat_B_26_min = rep(0, 50769),
                        climat_B_26_max = rep(0, 50769),
                        climat_B_85_mean = rep(0, 50769),
                        climat_B_85_min = rep(0, 50769),
                        climat_B_85_max = rep(0, 50769),
                        climat_C_26_mean = rep(0, 50769),
                        climat_C_26_min = rep(0, 50769),
                        climat_C_26_max = rep(0, 50769),
                        climat_C_85_mean = rep(0, 50769),
                        climat_C_85_min = rep(0, 50769),
                        climat_C_85_max = rep(0, 50769),
                        climat_D_26_mean = rep(0, 50769),
                        climat_D_26_min = rep(0, 50769),
                        climat_D_26_max = rep(0, 50769),
                        climat_D_85_mean = rep(0, 50769),
                        climat_D_85_min = rep(0, 50769),
                        climat_D_85_max = rep(0, 50769),
                        climat_E_26_mean = rep(0, 50769),
                        climat_E_26_min = rep(0, 50769),
                        climat_E_26_max = rep(0, 50769),
                        climat_E_85_mean = rep(0, 50769),
                        climat_E_85_min = rep(0, 50769),
                        climat_E_85_max = rep(0, 50769))



# Distribution by oceans
subset_Arc = subset(row.names(climats), climats$OCEAN == 'Arctic Ocean')
subset_Atl = subset(row.names(climats), climats$OCEAN == 'Atlantic Ocean') 
subset_Pac = subset(row.names(climats), climats$OCEAN == 'Pacific Ocean') 
subset_Ind = subset(row.names(climats), climats$OCEAN == 'Indian Ocean') 

oceans_Tg = data.frame(date= seq(as.Date("1961/1/1"), as.Date("2099/12/31"), "days"),
                        Arc_26_mean = rep(0, 50769),
                      Arc_26_min = rep(0, 50769),
                      Arc_26_max = rep(0, 50769),
                      Arc_85_mean = rep(0, 50769),
                      Arc_85_min = rep(0, 50769),
                      Arc_85_max = rep(0, 50769),
                        Atl_26_mean = rep(0, 50769),
                      Atl_26_min = rep(0, 50769),
                      Atl_26_max = rep(0, 50769),
                      Atl_85_mean = rep(0, 50769),
                      Atl_85_min = rep(0, 50769),
                      Atl_85_max = rep(0, 50769),
                        Pac_26_mean = rep(0, 50769),
                      Pac_26_min = rep(0, 50769),
                      Pac_26_max = rep(0, 50769),
                      Pac_85_mean = rep(0, 50769),
                      Pac_85_min = rep(0, 50769),
                      Pac_85_max = rep(0, 50769),
                        Ind_26_mean = rep(0, 50769),
                      Ind_26_min = rep(0, 50769),
                      Ind_26_max = rep(0, 50769),
                      Ind_85_mean = rep(0, 50769),
                      Ind_85_min = rep(0, 50769),
                      Ind_85_max = rep(0, 50769))


# Distribution by continent
Africa = subset(row.names(climats), climats$CONTINENT == 'Africa' & is.na(climats$OCEAN) == F)
Asia = subset(row.names(climats), climats$CONTINENT == 'Asia' & climats$OCEAN != 'Caspian Sea' & is.na(climats$OCEAN) == F) 
Europe = subset(row.names(climats), climats$CONTINENT == 'Europe' & climats$OCEAN != 'Caspian Sea' & is.na(climats$OCEAN) == F) 
N_Am = subset(row.names(climats), climats$CONTINENT == 'North America') 
S_Am = subset(row.names(climats), climats$CONTINENT == 'South America') 
Oceania = subset(row.names(climats), climats$CONTINENT == 'Oceania') 

continents_Tg = data.frame(date= seq(as.Date("1961/1/1"), as.Date("2099/12/31"), "days"),
                           Africa_26_mean = rep(0, 50769),
                           Africa_26_min = rep(0, 50769),
                           Africa_26_max = rep(0, 50769),
                           Africa_85_mean = rep(0, 50769),
                           Africa_85_min = rep(0, 50769),
                           Africa_85_max = rep(0, 50769),
                           Asia_26_mean = rep(0, 50769),
                           Asia_26_min = rep(0, 50769),
                           Asia_26_max = rep(0, 50769),
                           Asia_85_mean = rep(0, 50769),
                           Asia_85_min = rep(0, 50769),
                           Asia_85_max = rep(0, 50769),
                           Europe_26_mean = rep(0, 50769),
                           Europe_26_min = rep(0, 50769),
                           Europe_26_max = rep(0, 50769),
                           Europe_85_mean = rep(0, 50769),
                           Europe_85_min = rep(0, 50769),
                           Europe_85_max = rep(0, 50769),
                           N_Am_26_mean = rep(0, 50769),
                           N_Am_26_min = rep(0, 50769),
                           N_Am_26_max = rep(0, 50769),
                           N_Am_85_mean = rep(0, 50769),
                           N_Am_85_min = rep(0, 50769),
                           N_Am_85_max = rep(0, 50769),
                           S_Am_26_mean = rep(0, 50769),
                           S_Am_26_min = rep(0, 50769),
                           S_Am_26_max = rep(0, 50769),
                           S_Am_85_mean = rep(0, 50769),
                           S_Am_85_min = rep(0, 50769),
                           S_Am_85_max = rep(0, 50769),
                           Oceania_26_mean = rep(0, 50769),
                           Oceania_26_min = rep(0, 50769),
                           Oceania_26_max = rep(0, 50769),
                           Oceania_85_mean = rep(0, 50769),
                           Oceania_85_min = rep(0, 50769),
                           Oceania_85_max = rep(0, 50769))

for (file in files_climat){
  table = read.csv(file.path(path_final,file), header=T,sep=",")
  name = strsplit(file,".csv")
  print(name)
  
  name_short = str_sub(name,1,nchar(name)-4)
  
  table$date <- as.Date(with(table, paste(Y, M, D,sep="-")), "%Y-%m-%d")
  table <- table[order(table$date),]

  global_Tg$global_Tg_26_mean =   global_Tg$global_Tg_26_mean + table$DOC_flux_Tg_26_mean
  global_Tg$global_Tg_26_min = global_Tg$global_Tg_26_min + table$DOC_flux_Tg_26_min
  global_Tg$global_Tg_26_max = global_Tg$global_Tg_26_max + table$DOC_flux_Tg_26_max
  
  global_Tg$global_Tg_85_mean = global_Tg$global_Tg_85_mean + table$DOC_flux_Tg_85_mean
  global_Tg$global_Tg_85_min = global_Tg$global_Tg_85_min + table$DOC_flux_Tg_85_min
  global_Tg$global_Tg_85_max = global_Tg$global_Tg_85_max + table$DOC_flux_Tg_85_max

# Climates
      if (name %in% subset_A){
        n_A = n_A + 1
        climats_Tg$climat_A_26_mean = climats_Tg$climat_A_26_mean + table$DOC_flux_Tg_26_mean
        climats_Tg$climat_A_26_min = climats_Tg$climat_A_26_min + table$DOC_flux_Tg_26_min
        climats_Tg$climat_A_26_max = climats_Tg$climat_A_26_max + table$DOC_flux_Tg_26_max
        
        climats_Tg$climat_A_85_mean = climats_Tg$climat_A_85_mean + table$DOC_flux_Tg_85_mean
        climats_Tg$climat_A_85_min = climats_Tg$climat_A_85_min + table$DOC_flux_Tg_85_min
        climats_Tg$climat_A_85_max = climats_Tg$climat_A_85_max + table$DOC_flux_Tg_85_max
      } else if (name %in% subset_B) {
        n_B = n_B + 1
        climats_Tg$climat_B_26_mean = climats_Tg$climat_B_26_mean + table$DOC_flux_Tg_26_mean
        climats_Tg$climat_B_26_min = climats_Tg$climat_B_26_min + table$DOC_flux_Tg_26_min
        climats_Tg$climat_B_26_max = climats_Tg$climat_B_26_max + table$DOC_flux_Tg_26_max
        
        climats_Tg$climat_B_85_mean = climats_Tg$climat_B_85_mean + table$DOC_flux_Tg_85_mean
        climats_Tg$climat_B_85_min = climats_Tg$climat_B_85_min + table$DOC_flux_Tg_85_min
        climats_Tg$climat_B_85_max = climats_Tg$climat_B_85_max + table$DOC_flux_Tg_85_max        
      } else if (name %in% subset_C) {
        n_C = n_C + 1
        climats_Tg$climat_C_26_mean = climats_Tg$climat_C_26_mean + table$DOC_flux_Tg_26_mean
        climats_Tg$climat_C_26_min = climats_Tg$climat_C_26_min + table$DOC_flux_Tg_26_min
        climats_Tg$climat_C_26_max = climats_Tg$climat_C_26_max + table$DOC_flux_Tg_26_max
        
        climats_Tg$climat_C_85_mean = climats_Tg$climat_C_85_mean + table$DOC_flux_Tg_85_mean
        climats_Tg$climat_C_85_min = climats_Tg$climat_C_85_min + table$DOC_flux_Tg_85_min
        climats_Tg$climat_C_85_max = climats_Tg$climat_C_85_max + table$DOC_flux_Tg_85_max          
      } else if (name %in% subset_D) {
        n_D = n_D + 1
        climats_Tg$climat_D_26_mean = climats_Tg$climat_D_26_mean + table$DOC_flux_Tg_26_mean
        climats_Tg$climat_D_26_min = climats_Tg$climat_D_26_min + table$DOC_flux_Tg_26_min
        climats_Tg$climat_D_26_max = climats_Tg$climat_D_26_max + table$DOC_flux_Tg_26_max
        
        climats_Tg$climat_D_85_mean = climats_Tg$climat_D_85_mean + table$DOC_flux_Tg_85_mean
        climats_Tg$climat_D_85_min = climats_Tg$climat_D_85_min + table$DOC_flux_Tg_85_min
        climats_Tg$climat_D_85_max = climats_Tg$climat_D_85_max + table$DOC_flux_Tg_85_max         
      } else if (name %in% subset_E) {
        n_E = n_E + 1
        climats_Tg$climat_E_26_mean = climats_Tg$climat_E_26_mean + table$DOC_flux_Tg_26_mean
        climats_Tg$climat_E_26_min = climats_Tg$climat_E_26_min + table$DOC_flux_Tg_26_min
        climats_Tg$climat_E_26_max = climats_Tg$climat_E_26_max + table$DOC_flux_Tg_26_max
        
        climats_Tg$climat_E_85_mean = climats_Tg$climat_E_85_mean + table$DOC_flux_Tg_85_mean
        climats_Tg$climat_E_85_min = climats_Tg$climat_E_85_min + table$DOC_flux_Tg_85_min
        climats_Tg$climat_E_85_max = climats_Tg$climat_E_85_max + table$DOC_flux_Tg_85_max         
      }
        
# Oceans        
      if (name %in% subset_Arc){
        oceans_Tg$Arc_26_mean = oceans_Tg$Arc_26_mean + table$DOC_flux_Tg_26_mean
        oceans_Tg$Arc_26_min = oceans_Tg$Arc_26_min + table$DOC_flux_Tg_26_min
        oceans_Tg$Arc_26_max = oceans_Tg$Arc_26_max + table$DOC_flux_Tg_26_max
        
        oceans_Tg$Arc_85_mean = oceans_Tg$Arc_85_mean + table$DOC_flux_Tg_85_mean
        oceans_Tg$Arc_85_min = oceans_Tg$Arc_85_min + table$DOC_flux_Tg_85_min
        oceans_Tg$Arc_85_max = oceans_Tg$Arc_85_max + table$DOC_flux_Tg_85_max
      } else if (name %in% subset_Atl) {
        oceans_Tg$Atl_26_mean = oceans_Tg$Atl_26_mean + table$DOC_flux_Tg_26_mean
        oceans_Tg$Atl_26_min = oceans_Tg$Atl_26_min + table$DOC_flux_Tg_26_min
        oceans_Tg$Atl_26_max = oceans_Tg$Atl_26_max + table$DOC_flux_Tg_26_max
        
        oceans_Tg$Atl_85_mean = oceans_Tg$Atl_85_mean + table$DOC_flux_Tg_85_mean
        oceans_Tg$Atl_85_min = oceans_Tg$Atl_85_min + table$DOC_flux_Tg_85_min
        oceans_Tg$Atl_85_max = oceans_Tg$Atl_85_max + table$DOC_flux_Tg_85_max        
      } else if (name %in% subset_Pac) {
        oceans_Tg$Pac_26_mean = oceans_Tg$Pac_26_mean + table$DOC_flux_Tg_26_mean
        oceans_Tg$Pac_26_min = oceans_Tg$Pac_26_min + table$DOC_flux_Tg_26_min
        oceans_Tg$Pac_26_max = oceans_Tg$Pac_26_max + table$DOC_flux_Tg_26_max
        
        oceans_Tg$Pac_85_mean = oceans_Tg$Pac_85_mean + table$DOC_flux_Tg_85_mean
        oceans_Tg$Pac_85_min = oceans_Tg$Pac_85_min + table$DOC_flux_Tg_85_min
        oceans_Tg$Pac_85_max = oceans_Tg$Pac_85_max + table$DOC_flux_Tg_85_max          
      } else if (name %in% subset_Ind) {
        oceans_Tg$Ind_26_mean = oceans_Tg$Ind_26_mean + table$DOC_flux_Tg_26_mean
        oceans_Tg$Ind_26_min = oceans_Tg$Ind_26_min + table$DOC_flux_Tg_26_min
        oceans_Tg$Ind_26_max = oceans_Tg$Ind_26_max + table$DOC_flux_Tg_26_max
        
        oceans_Tg$Ind_85_mean = oceans_Tg$Ind_85_mean + table$DOC_flux_Tg_85_mean
        oceans_Tg$Ind_85_min = oceans_Tg$Ind_85_min + table$DOC_flux_Tg_85_min
        oceans_Tg$Ind_85_max = oceans_Tg$Ind_85_max + table$DOC_flux_Tg_85_max
      }        
 
# Continents
      if (name %in% Africa){
        continents_Tg$Africa_26_mean = continents_Tg$Africa_26_mean + table$DOC_flux_Tg_26_mean
        continents_Tg$Africa_26_min = continents_Tg$Africa_26_min + table$DOC_flux_Tg_26_min
        continents_Tg$Africa_26_max = continents_Tg$Africa_26_max + table$DOC_flux_Tg_26_max
        
        continents_Tg$Africa_85_mean = continents_Tg$Africa_85_mean + table$DOC_flux_Tg_85_mean
        continents_Tg$Africa_85_min = continents_Tg$Africa_85_min + table$DOC_flux_Tg_85_min
        continents_Tg$Africa_85_max = continents_Tg$Africa_85_max + table$DOC_flux_Tg_85_max
      } else if (name %in% Asia) {
        continents_Tg$Asia_26_mean = continents_Tg$Asia_26_mean + table$DOC_flux_Tg_26_mean
        continents_Tg$Asia_26_min = continents_Tg$Asia_26_min + table$DOC_flux_Tg_26_min
        continents_Tg$Asia_26_max = continents_Tg$Asia_26_max + table$DOC_flux_Tg_26_max
        
        continents_Tg$Asia_85_mean = continents_Tg$Asia_85_mean + table$DOC_flux_Tg_85_mean
        continents_Tg$Asia_85_min = continents_Tg$Asia_85_min + table$DOC_flux_Tg_85_min
        continents_Tg$Asia_85_max = continents_Tg$Asia_85_max + table$DOC_flux_Tg_85_max        
      } else if (name %in% Europe) {
        continents_Tg$Europe_26_mean = continents_Tg$Europe_26_mean + table$DOC_flux_Tg_26_mean
        continents_Tg$Europe_26_min = continents_Tg$Europe_26_min + table$DOC_flux_Tg_26_min
        continents_Tg$Europe_26_max = continents_Tg$Europe_26_max + table$DOC_flux_Tg_26_max
        
        continents_Tg$Europe_85_mean = continents_Tg$Europe_85_mean + table$DOC_flux_Tg_85_mean
        continents_Tg$Europe_85_min = continents_Tg$Europe_85_min + table$DOC_flux_Tg_85_min
        continents_Tg$Europe_85_max = continents_Tg$Europe_85_max + table$DOC_flux_Tg_85_max          
      } else if (name %in% N_Am) {
        continents_Tg$N_Am_26_mean = continents_Tg$N_Am_26_mean + table$DOC_flux_Tg_26_mean
        continents_Tg$N_Am_26_min = continents_Tg$N_Am_26_min + table$DOC_flux_Tg_26_min
        continents_Tg$N_Am_26_max = continents_Tg$N_Am_26_max + table$DOC_flux_Tg_26_max
        
        continents_Tg$N_Am_85_mean = continents_Tg$N_Am_85_mean + table$DOC_flux_Tg_85_mean
        continents_Tg$N_Am_85_min = continents_Tg$N_Am_85_min + table$DOC_flux_Tg_85_min
        continents_Tg$N_Am_85_max = continents_Tg$N_Am_85_max + table$DOC_flux_Tg_85_max
      } else if (name %in% S_Am) {
        continents_Tg$S_Am_26_mean = continents_Tg$S_Am_26_mean + table$DOC_flux_Tg_26_mean
        continents_Tg$S_Am_26_min = continents_Tg$S_Am_26_min + table$DOC_flux_Tg_26_min
        continents_Tg$S_Am_26_max = continents_Tg$S_Am_26_max + table$DOC_flux_Tg_26_max
        
        continents_Tg$S_Am_85_mean = continents_Tg$S_Am_85_mean + table$DOC_flux_Tg_85_mean
        continents_Tg$S_Am_85_min = continents_Tg$S_Am_85_min + table$DOC_flux_Tg_85_min
        continents_Tg$S_Am_85_max = continents_Tg$S_Am_85_max + table$DOC_flux_Tg_85_max
      } else if (name_short %in% Oceania) {
        continents_Tg$Oceania_26_mean = continents_Tg$Oceania_26_mean + table$DOC_flux_Tg_26_mean
        continents_Tg$Oceania_26_min = continents_Tg$Oceania_26_min + table$DOC_flux_Tg_26_min
        continents_Tg$Oceania_26_max = continents_Tg$Oceania_26_max + table$DOC_flux_Tg_26_max
        
        continents_Tg$Oceania_85_mean = continents_Tg$Oceania_85_mean + table$DOC_flux_Tg_85_mean
        continents_Tg$Oceania_85_min = continents_Tg$Oceania_85_min + table$DOC_flux_Tg_85_min
        continents_Tg$Oceania_85_max = continents_Tg$Oceania_85_max + table$DOC_flux_Tg_85_max
      }    
}
n_tot = n_A + n_B + n_C + n_D + n_E

write.csv(global_Tg, file="global_Tg.csv")
write.csv(climats_Tg, file="climats_Tg.csv")
write.csv(oceans_Tg, file="oceans_Tg.csv")
write.csv(continents_Tg, file="continents_Tg.csv")


global_Tg$julian <- as.POSIXlt(global_Tg$date)$yday + 1
global_Tg_past <- global_Tg[year(global_Tg$date) >= 1971 & year(global_Tg$date) <= 2000,]
results_global_Tg_past = as.data.frame(global_Tg_past[,2:8] %>%
                                   group_by(julian) %>%
                                   summarise_each(funs(mean)))[,2:7]

global_Tg_fut <- global_Tg[year(global_Tg$date) >= 2071,]
results_global_Tg_fut = as.data.frame(global_Tg_fut[,2:8] %>%
                                                 group_by(julian) %>%
                                                 summarise_each(funs(mean)))[,2:7]

climats_Tg$julian <- as.POSIXlt(climats_Tg$date)$yday + 1
climats_Tg_past <- climats_Tg[year(climats_Tg$date) >= 1971 & year(climats_Tg$date) <= 2000,]
results_climats_Tg_past = as.data.frame(climats_Tg_past[,2:32] %>%
                                                 group_by(julian) %>%
                                                 summarise_each(funs(mean)))[,2:31]
climats_Tg_fut <- climats_Tg[year(climats_Tg$date) >= 2071,]
results_climats_Tg_fut = as.data.frame(climats_Tg_fut[,2:32] %>%
                                                group_by(julian) %>%
                                                summarise_each(funs(mean)))[,2:31]

oceans_Tg$julian <- as.POSIXlt(oceans_Tg$date)$yday + 1
oceans_Tg_past <- oceans_Tg[year(oceans_Tg$date) >= 1971 & year(oceans_Tg$date) <= 2000,]
results_oceans_Tg_past = as.data.frame(oceans_Tg_past[,2:26] %>%
                                                  group_by(julian) %>%
                                                  summarise_each(funs(mean)))[,2:25]
colnames(results_oceans_Tg_past)=c("sum26_2","sum26_2_min","sum26_2_max","sum85_2","sum85_2_min","sum85_2_max",
                                  "sum26_1","sum26_1_min","sum26_1_max","sum85_1","sum85_1_min","sum85_1_max",
                                  "sum26_3","sum26_3_min","sum26_3_max","sum85_3","sum85_3_min","sum85_3_max",
                                  "sum26_4","sum26_4_min","sum26_4_max","sum85_4","sum85_4_min","sum85_4_max")

oceans_Tg_fut <- oceans_Tg[year(oceans_Tg$date) >= 2071,]
results_oceans_Tg_fut = as.data.frame(oceans_Tg_fut[,2:26] %>%
                                                  group_by(julian) %>%
                                                  summarise_each(funs(mean)))[,2:25]
colnames(results_oceans_Tg_fut)=c("sum26_2","sum26_2_min","sum26_2_max","sum85_2","sum85_2_min","sum85_2_max",
                                  "sum26_1","sum26_1_min","sum26_1_max","sum85_1","sum85_1_min","sum85_1_max",
                                  "sum26_3","sum26_3_min","sum26_3_max","sum85_3","sum85_3_min","sum85_3_max",
                                  "sum26_4","sum26_4_min","sum26_4_max","sum85_4","sum85_4_min","sum85_4_max")

continents_Tg$julian <- as.POSIXlt(continents_Tg$date)$yday + 1
continents_Tg_past <- continents_Tg[year(continents_Tg$date) >= 1971 & year(continents_Tg$date) <= 2000,]
results_continents_Tg_past = as.data.frame(continents_Tg_past[,2:38] %>%
                                                    group_by(julian) %>%
                                                    summarise_each(funs(mean)))[,2:37]
continents_Tg_fut <- continents_Tg[year(continents_Tg$date) >= 2071,]
results_continents_Tg_fut = as.data.frame(continents_Tg_fut[,2:38] %>%
                                                 group_by(julian) %>%
                                                 summarise_each(funs(mean)))[,2:37]


# Figure 5 ####
png(filename=paste("Flux_oceans_26.png",sep=""),
    width = 2800, height = 2800,  res = 300)
ggplot(data = results_oceans_Tg_fut, aes(x = 1:nrow(results_oceans_Tg_fut))) +
  xlab("") +
  ylab("Average daily DOC export (TgC."~day^-1*")") +
  expand_limits(x=0, y=0) +
  geom_line(aes(y=sum26_1, color='sum26_1'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum26_1_min, ymax=sum26_1_max, fill="sum26_1", alpha="sum26_1"), show.legend = FALSE) +
  geom_line(aes(y=sum26_2, color='sum26_2'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum26_2_min, ymax=sum26_2_max, fill="sum26_2", alpha="sum26_2"), show.legend = FALSE) +
  geom_line(aes(y=sum26_3, color='sum26_3'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum26_3_min, ymax=sum26_3_max, fill="sum26_3", alpha="sum26_3"), show.legend = FALSE) +
  geom_line(aes(y=sum26_4, color='sum26_4'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum26_4_min, ymax=sum26_4_max, fill="sum26_4", alpha="sum26_4"), show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 28, angle = 90,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size = 28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 20),
        axis.line = element_line(colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size=1),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_linetype_manual(name="Ocean supplied", values = c(rep("solid", 1), rep("solid", 1),
                                                          rep("dotted", 1), rep("dashed", 1)),
                        labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_color_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_fill_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                    labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  theme(legend.justification=c(0,1), legend.position=c(0.7,0.96),
        legend.background = element_rect(fill="white", size=.5),
        legend.title=element_text(size=19), legend.text=element_text(size = 19),
        legend.key.height=unit(2,"line"), legend.key.width=unit(2,"line"))+
  scale_x_continuous(limits = c(0,366), breaks = c(15,45,76,106,136,167,197,228,259,289,320,350),
                     minor_breaks = c(90,181,273), labels = c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0, 0.5, 0.1),
                     minor_breaks = seq(0, 0.5, 0.1))
dev.off()

png(filename=paste("Flux_oceans_85.png",sep=""),
    width = 2800, height = 2800,  res = 300)
ggplot(data = results_oceans_Tg_fut, aes(x = 1:nrow(results_oceans_Tg_fut))) +
  xlab("") +
  ylab("Average daily DOC export (TgC."~day^-1*")") +
  expand_limits(x=0, y=0) +
  geom_line(aes(y=sum85_1, color='sum85_1'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum85_1_min, ymax=sum85_1_max, fill="std_1", alpha="std_1"), show.legend = FALSE) +
  geom_line(aes(y=sum85_2, color='sum85_2'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum85_2_min, ymax=sum85_2_max, fill="std_2", alpha="std_2"), show.legend = FALSE) +
  geom_line(aes(y=sum85_3, color='sum85_3'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum85_3_min, ymax=sum85_3_max, fill="std_3", alpha="std_3"), show.legend = FALSE) +
  geom_line(aes(y=sum85_4, color='sum85_4'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(results_oceans_Tg_fut), ymin=sum85_4_min, ymax=sum85_4_max, fill="std_4", alpha="std_4"), show.legend = FALSE) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 28, angle = 90,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size = 28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 20),
        axis.line = element_line(colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size=1),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_linetype_manual(name="Ocean supplied", values = c(rep("solid", 1), rep("solid", 1),
                                                          rep("dotted", 1), rep("dashed", 1)),
                        labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_color_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_fill_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                    labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  theme(legend.justification=c(0,1), legend.position=c(0.7,0.96),
        legend.background = element_rect(fill="white", size=.5),
        legend.title=element_text(size=19), legend.text=element_text(size = 19),
        legend.key.height=unit(2,"line"), legend.key.width=unit(2,"line"))+
  scale_x_continuous(limits = c(0,366), breaks = c(15,45,76,106,136,167,197,228,259,289,320,350),
                     minor_breaks = c(90,181,273), labels = c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0, 0.5, 0.1),
                     minor_breaks = seq(0, 0.5, 0.1))
dev.off()


# Anomalies
ano_oceans <- 100*(results_oceans_Tg_fut-results_oceans_Tg_past)/results_oceans_Tg_past

png(filename="Ano_oceans_26.png",
    width = 2800, height = 2000,  res = 300)
ggplot(data = ano_oceans, aes(x = 1:nrow(ano_oceans))) +
  xlab("") +
  ylab(expression(~F[DOC] * " change (%)")) +
  expand_limits(x=0, y=0) +
  geom_line(aes(y=sum26_1, color='sum26_1'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum26_1_min, ymax=sum26_1_max, fill="std_1", alpha="std_1"), show.legend = FALSE) +
  geom_line(aes(y=sum26_2, color='sum26_2'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum26_2_min, ymax=sum26_2_max, fill="std_2", alpha="std_2"), show.legend = FALSE) +
  geom_line(aes(y=sum26_3, color='sum26_3'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum26_3_min, ymax=sum26_3_max, fill="std_3", alpha="std_3"), show.legend = FALSE) +
  geom_line(aes(y=sum26_4, color='sum26_4'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum26_4_min, ymax=sum26_4_max, fill="std_4", alpha="std_4"), show.legend = FALSE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed" , size=1.2) +
  theme(axis.title.y = element_text(size = 28, angle = 90,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size = 28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 20),
        axis.line = element_line(colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size=1),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_linetype_manual(name="Ocean supplied", values = c(rep("solid", 1), rep("solid", 1),
                                                          rep("dotted", 1), rep("dashed", 1)),
                        labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_color_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_fill_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                    labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  theme(legend.justification=c(0,1), legend.position="none",
        legend.background = element_rect(fill="white", size=.5),
        legend.title=element_text(size=19), legend.text=element_text(size = 19),
        legend.key.height=unit(2,"line"), legend.key.width=unit(2,"line"))+
  scale_x_continuous(limits = c(0,366), breaks = c(15,45,76,106,136,167,197,228,259,289,320,350),
                     minor_breaks = c(90,181,273), labels = c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  scale_y_continuous(limits = c(-100,700), breaks = seq(-100,700, 100),
                     minor_breaks = seq(-100,150, 50))
dev.off()

png(filename="E:/Papiers/Papier Global carbon futur/Ano_oceans_85.png",
    width = 2800, height = 2000,  res = 300)
ggplot(data = ano_oceans, aes(x = 1:nrow(ano_oceans))) +
  xlab("") +
  ylab(expression(~F[DOC] * " change (%)")) +
  expand_limits(x=0, y=0) +
  geom_line(aes(y=sum85_1, color='sum85_1'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum85_1_min, ymax=sum85_1_max, fill="std_1", alpha="std_1"), show.legend = FALSE) +
  geom_line(aes(y=sum85_2, color='sum85_2'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum85_2_min, ymax=sum85_2_max, fill="std_2", alpha="std_2"), show.legend = FALSE) +
  geom_line(aes(y=sum85_3, color='sum85_3'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum85_3_min, ymax=sum85_3_max, fill="std_3", alpha="std_3"), show.legend = FALSE) +
  geom_line(aes(y=sum85_4, color='sum85_4'), size=1.2) +
  geom_ribbon(aes(x=1:nrow(ano_oceans), ymin=sum85_4_min, ymax=sum85_4_max, fill="std_4", alpha="std_4"), show.legend = FALSE) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed" , size=1.2) +
  theme(axis.title.y = element_text(size = 28, angle = 90,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(size = 28, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.y = element_text(size = 25), axis.text.x = element_text(size = 20),
        axis.line = element_line(colour = "black", size=1),
        axis.ticks = element_line(colour = "black", size=1),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank()) + #element_line(colour = "darkgrey", size=0.5)) +
  scale_linetype_manual(name="Ocean supplied", values = c(rep("solid", 1), rep("solid", 1),
                                                          rep("dotted", 1), rep("dashed", 1)),
                        labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_color_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_fill_manual(values=c('darkgreen','blue','darkorange4','goldenrod1'), name="Ocean supplied",
                    labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  scale_alpha_manual(values=c(0.3,0.3,0.3,0.3), name="Ocean supplied",
                     labels=c("Atlantic Ocean","Arctic Ocean","Pacific Ocean","Indian Ocean")) +
  theme(legend.justification=c(0,1), legend.position="none",
        legend.background = element_rect(fill="white", size=.5),
        legend.title=element_text(size=19), legend.text=element_text(size = 19),
        legend.key.height=unit(2,"line"), legend.key.width=unit(2,"line"))+
  scale_x_continuous(limits = c(0,366), breaks = c(15,45,76,106,136,167,197,228,259,289,320,350),
                     minor_breaks = c(90,181,273), labels = c("J","F","M","A","M","J","J","A","S","O","N","D")) +
  scale_y_continuous(limits = c(-100,700), breaks = seq(-100,700, 100),
                     minor_breaks = seq(-100,150, 50))
dev.off()


# Figure 6 ####
basins=read_sf(layer = "basins_to_oceans")
basins <- st_zm(basins, drop=T)

countries=read_sf(dlayer="Countries_without_antartica")
countries <- st_zm(countries, drop=T)


tot_26 = read.csv("flux_DOC_by_rivers_Tg_26.csv", sep=";")
tot_85 = read.csv("flux_DOC_by_rivers_Tg_85.csv")

map_countries = tm_shape(countries) + tm_fill()

diff_26 = 100* (tot_26[,14:23] - tot_26[,4:13])/tot_26[,4:13]
diff_26$river = tot_26$river

basins_26=merge(basins, diff_26, by.x = "wribasin_u", by.y = "river")

map_26 = map_countries + tm_shape(basins_26) + tm_borders() +
  tm_fill(col="total_mean_fut", breaks=c(-100,-10,-5,0,5,10,100),
          labels=c("< -10","-10 to -5","-5 to 0",
                   "0 to 5","5 to 10","> 10"),
          title=expression(~F[DOC] * " change (%)"),
          palette="RdYlGn") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_26, dpi=300,filename = "change_map_26.png")

diff_85 = 100* (tot_85[,14:23] - tot_85[,4:13])/tot_85[,4:13]
diff_85$river = tot_85$river

basins_85=merge(basins, diff_85, by.x = "wribasin_u", by.y = "river")
map_85 = map_countries + tm_shape(basins_85) + tm_borders() +
  tm_fill(col="total_mean_fut", breaks=c(-100,-10,-5,0,5,10,100),
          labels=c("< -10","-10 to -5","-5 to 0",
                   "0 to 5","5 to 10","> 10"),
          title=expression(~F[DOC] * " change (%)"),
          palette="RdYlGn") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_85, dpi=300,filename = "change_map_85.png")



### Map total
basins_26_past=merge(basins, tot_26, by.x = "wribasin_u", by.y = "river")

map_26_past = map_countries + tm_shape(basins_26_past) + tm_borders() +
  tm_fill(col="total_mean_past", breaks=c(0,0.5,1.5,3.5,16,30),
          title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
          palette="YlOrRd") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_26_past, dpi=300,filename = "map_26_past.png")

basins_85_past=merge(basins, tot_85, by.x = "wribasin_u", by.y = "river")

map_85_past = map_countries + tm_shape(basins_85_past) + tm_borders() +
  tm_fill(col="total_mean_past", breaks=c(0,0.5,1.5,3.5,16,30),
          title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
          palette="YlOrRd") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_85_past, dpi=300,filename = "map_85_past.png")



basins_26_fut=merge(basins, tot_26, by.x = "wribasin_u", by.y = "river")

map_26_fut = map_countries + tm_shape(basins_26_fut) + tm_borders() +
  tm_fill(col="total_mean_fut", breaks=c(0,0.5,1.5,3.5,16,30),
          title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
          palette="YlOrRd") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_26_fut, dpi=300,filename = "map_26_fut.png")

basins_85_fut=merge(basins, tot_85, by.x = "wribasin_u", by.y = "river")

map_85_fut = map_countries + tm_shape(basins_85_fut) + tm_borders() +
  tm_fill(col="total_mean_fut", breaks=c(0,0.5,1.5,3.5,16,30),
          title=expression(~F[DOC] * " (TgC."~yr^-1*")"),
          palette="YlOrRd") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_85_fut, dpi=300,filename = "map_85_fut.png")


# Appendix 1 & 2 ####
climats_map = climats
row.names(climats_map) = seq(1,nrow(climats_map),1)
climats_map$river = row.names(climats)
basins_climats <- merge(basins, climats_map, by.x = "wribasin_u", by.y = "river")

basins_climats$change_soc_26 = 100*(basins_climats$soc_26_70 - basins_climats$soc_past)/basins_climats$soc_past
basins_climats$change_soc_85 = 100*(basins_climats$soc_85_70 - basins_climats$soc_past)/basins_climats$soc_past

basins_climats$change_T_26 = (basins_climats$T_26_70 - basins_climats$T_past)
basins_climats$change_T_85 = (basins_climats$T_85_70 - basins_climats$T_past)

basins_climats$P_past_year = 12*basins_climats$P_past
basins_climats$P_26_70_year = 12*basins_climats$P_26_70
basins_climats$P_85_70_year = 12*basins_climats$P_85_70

basins_climats$change_P_26 = (basins_climats$P_26_70_year - basins_climats$P_past_year)
basins_climats$change_P_85 = (basins_climats$P_85_70_year - basins_climats$P_past_year)


# Discharge
map_Q_26 = map_countries + tm_shape(basins_26) + tm_borders() +
  tm_fill(col="total_Q_fut", breaks=c(-100,-10,-5,0,5,10,100),
          labels=c("< -10","-10 to -5","-5 to 0",
                   "0 to 5","5 to 10","> 10"),
          title=expression("Discharge change (%)"),
          palette="RdBu") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_Q_26, dpi=300,filename = "change_Q_map_26.png")

map_Q_85 = map_countries + tm_shape(basins_85) + tm_borders() +
  tm_fill(col="total_Q_fut", breaks=c(-100,-10,-5,0,5,10,100),
          labels=c("< -10","-10 to -5","-5 to 0",
                   "0 to 5","5 to 10","> 10"),
          title=expression("Discharge change (%)"),
          palette="RdBu") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_Q_85, dpi=300,filename = "change_Q_map_85.png")


basins_26_tot=merge(basins, tot_26, by.x = "wribasin_u", by.y = "river")
map_Qtot_past = map_countries + tm_shape(basins_26_tot) + tm_borders() +
  tm_fill(col="total_Q_past", breaks=c(0,100,200,500,1000,10000),
          labels=c("< 100", "100 to 200", "200 to 500", "500 to 1000", "> 1000"),
          title=expression("Discharge ("~km^3*""~yr^-1*")"),
          palette="Blues") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_Qtot_past, dpi=300,filename = "tot_Q_map_past.png")

# Temperature
map_T_past = map_countries + tm_shape(basins_climats) + tm_borders() +
  tm_fill(col="T_past", breaks=c(-20,-10,-5,0,5,10,15,20,25,30),
          labels=c("< -10", "-10 to -5", "-5 to 0", "0 to 5", "5 to 10", "10 to 15",
                   "15 to 20", "20 to 25", "25 to 30"),
          title=expression("Air temperature (°C)"),
          palette="-RdYlBu") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_T_past, dpi=300,filename = "T_map_past.png")

map_T_26 = map_countries + tm_shape(basins_climats) + tm_borders() +
  tm_fill(col="change_T_26", breaks=c(-10,0,1,2,5,10),
          labels=c("-10 to 0", "0 to 1", "1 to 2", "2 to 5", "5 to 10"),
          title=expression("Temperature deviation (°C)"),
          palette="-RdYlBu") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_T_26, dpi=300,filename = "T_map_26.png")

map_T_85 = map_countries + tm_shape(basins_climats) + tm_borders() +
  tm_fill(col="change_T_85", breaks=c(-10,0,1,2,5,10),
          labels=c("-10 to 0", "0 to 1", "1 to 2", "2 to 5", "5 to 10"),
          title=expression("Temperature deviation (°C)"),
          palette="-RdYlBu") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_T_85, dpi=300,filename = "T_map_85.png")

# Precipitation
map_P_past = map_countries + tm_shape(basins_climats) + tm_borders() +
  tm_fill(col="P_past_year", breaks=c(0,200,400,600,800,1000,10000),
          labels=c("0 to 200", "200 to 400", "400 to 600", "600 to 800", "800 to 1000", "> 1000"),
          title=expression("Precipitation (mm"~yr^-1*")"),
          palette="Blues") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.5,
            legend.text.size = 1.2)
tmap_save(tm = map_P_past, dpi=300,filename = "P_map_past.png")

map_P_26 = map_countries + tm_shape(basins_climats) + tm_borders() +
  tm_fill(col="change_P_26", breaks=c(-600,-200,-100,0,100,200,500,1000,10000),
          labels=c("< -200", "-200 to -100", "-100 to 0", "0 to 100", "100 to 200", "200 to 500", "500 to 1000", "> 1000"),
          title=expression("Precipitation deviation (mm"~yr^-1*")"),
          palette="RdBu") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.3,
            legend.text.size = 1.2)
tmap_save(tm = map_P_26, dpi=300,filename = "P_map_26.png")

map_P_85 = map_countries + tm_shape(basins_climats) + tm_borders() +
  tm_fill(col="change_P_85", breaks=c(-600,-200,-100,0,100,200,500,1000,10000),
          labels=c("< -200", "-200 to -100", "-100 to 0", "0 to 100", "100 to 200", "200 to 500", "500 to 1000", "> 1000"),
          title=expression("Precipitation deviation (mm"~yr^-1*")"),
          palette="RdBu") +
  tm_layout(frame=FALSE,legend.position = c("left","bottom"),
            legend.title.size = 1.3,
            legend.text.size = 1.2)
tmap_save(tm = map_P_85, dpi=300,filename = "P_map_85.png")


