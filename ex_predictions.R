# libraries
library(tidyr)
library(dplyr)
library(rstan)
library(ggplot2)
library(sf)
library(pROC)

# focal species = Fagus sylvatica
sp="Fagus sylvatica"

load(paste0("model_fit_safetymargin/",sp,"_2sm.RData"))
sp.trait=readRDS("data/df.species") |> 
  filter(species==sp)

occurence<-readRDS("data/occurence") |> 
  filter(species==sp) |> 
  mutate(hsm=(psi_eraday_real/1000)-sp.trait$px,
         fsm=tmin_era-sp.trait$lt50)

# get pars posteriors
post<-as.data.frame(fit.2var)

post |> 
  pivot_longer(cols=everything()) |> 
  ggplot(aes(value))+
  geom_density()+
  facet_wrap(~name,ncol=1,scales="free")
post |> 
  mutate(nsamp=row_number()) |> 
  pivot_longer(cols=-nsamp) |> 
  ggplot(aes(nsamp,value))+
  geom_line()+
  facet_wrap(~name,ncol=1,scales="free_y")

# plot map of Fagus presence 
worldmap <- st_as_sf(rworldmap::getMap(resolution = "high"))

occurence |>
  group_by(presence) |> 
  sample_n(10000) |> 
  ggplot()+
  geom_point(aes(x=x,y=y,color=as.factor(presence)))+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)

# plot safety margin
occurence |> 
  filter(!is.na(hsm),!is.na(fsm)) |> 
  mutate(hsm=cut(hsm,
                 breaks=c(-Inf, -15, -10, -5,-3, -2, -1, 0,1,2,3,5,Inf)),
         fsm=cut(fsm,
                  breaks=c(-Inf, -50,-40, -30,-20,-15,-10,-5,0,5,10,Inf))
         ) |> 
  sample_n(20000) |> 
  ggplot(aes(x=x,y=y,color=hsm))+
  geom_point()+
  scale_color_brewer(palette="RdYlBu")
occurence |> 
  mutate(hsm=(psi_eraday_real/1000)-sp.trait$px,
         fsm=tmin_era-sp.trait$lt50) |> 
  filter(!is.na(hsm),!is.na(fsm)) |> 
  mutate(hsm=cut(hsm,
                 breaks=c(-Inf, -15, -10, -5,-3, -2, -1, 0,1,2,3,5,Inf)),
         fsm=cut(fsm,
                 breaks=c(-Inf, -50,-40, -30,-20,-15,-10,-5,0,5,10,Inf))
  ) |> 
  sample_n(20000) |> 
  ggplot(aes(x=x,y=y,color=fsm))+
  geom_point()+
  scale_color_brewer(palette="RdYlBu")


# trajectory with uncertainty
fsm.05<- quantile(occurence$fsm,prob=0.05,na.rm=TRUE)[[1]]
fsm.95<- quantile(occurence$fsm,prob=0.95,na.rm=TRUE)[[1]]
hsm.05<- quantile(occurence$hsm,prob=0.05,na.rm=TRUE)[[1]]
hsm.95<- quantile(occurence$hsm,prob=0.95,na.rm=TRUE)[[1]]


post |> 
  crossing(data.frame(xsm_val=c(seq(fsm.05,
                                    fsm.95,
                                    length.out=100),
                                seq(hsm.05,
                                    hsm.95,
                                    length.out=100)),
                      xsm_name=c(rep("fsm",100),
                                 rep("hsm",100)),
                      xsm_name_l=c(rep("Frost safety margins (°C)",100),
                                   rep("Hydraulic safety margins (MPa)",100)))) |> 
  mutate(pred=case_when(xsm_name=="fsm"~K_int/((1+exp(-r_hsm*(hsm.95-t_hsm)))*
                                                  (1+exp(-r_fsm*(xsm_val-t_fsm)))),
                        xsm_name=="hsm"~K_int/((1+exp(-r_fsm*(fsm.95-t_fsm)))*
                                                 (1+exp(-r_hsm*(xsm_val-t_hsm)))))) |> 
  group_by(xsm_val,xsm_name_l) |> 
  summarise(med_pred=median(pred),
            q05_pred=quantile(pred,prob=0.05),
            q95_pred=quantile(pred,prob=0.95)) |> 
  ggplot()+
  geom_line(aes(xsm_val,med_pred))+
  geom_ribbon(aes(x=xsm_val,ymin=q05_pred,ymax=q95_pred),alpha=0.2)+
  facet_wrap(~xsm_name_l,scales = "free")+
  theme_minimal()+
  theme(axis.title= element_blank(),
        axis.ticks.x = element_line(linewidth = 0.6),
        axis.line = element_line(linewidth = 0.6),
        legend.position = "none",
        # panel.grid = element_blank(),
        # panel.grid.major = element_line(color = "lightgrey",
        #                                   size =0.02),
        axis.text.y = element_text(hjust=0.5),
        strip.placement = "outside",
        strip.text = element_text(size=11,vjust=-1.1),
        text=element_text(size=11))


# make predictions 

## determine the accurate threshold
post_med<-apply(post,MARGIN=2,median)
tss<-occurence |> 
  mutate(pred=post_med[["K_int"]]/
           ((1+exp(-post_med[["r_fsm"]]*(fsm-post_med[["t_fsm"]])))*
              (1+exp(-post_med[["r_hsm"]]*(hsm-post_med[["t_hsm"]]))))) |> 
  filter(!is.na(fsm),!is.na(hsm))


calc_tss <- function(threshold, observed, predicted_probs) {
  observed=as.factor(observed)
  predicted <- ifelse(predicted_probs > threshold, 1, 0)
  predicted=factor(predicted,levels=c(0,1))
  conf_matrix <- table(observed, predicted)
  TP <- conf_matrix[2, 2]
  FN <- conf_matrix[2, 1]
  TN <- conf_matrix[1, 1]
  FP <- conf_matrix[1, 2]
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  TSS <- sensitivity + specificity - 1
  return(TSS)
}
observed <- tss$presence
pred<- tss$pred
threshold_max<-post_med[["K_int"]]
thresholds <- seq(0,threshold_max, length.out=100)
tss_sfm <- sapply(thresholds, calc_tss, observed, pred)


opt_prob=thresholds[which.max(tss_sfm)]


## make prediction and plot
tss |> 
  mutate(pred_presence=pred>opt_prob) |> 
  group_by(presence) |> 
  sample_n(10000) |> 
  pivot_longer(cols=c("presence","pred_presence")) |> 
  ggplot()+
  geom_point(aes(x=x,y=y,color=value))+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)+
  facet_wrap(~name)
