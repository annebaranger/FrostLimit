##phenofit prediction

## 0. Setup
# Vérification que les packages nécessaires soient installés
list.of.packages <- c("data.table", "terra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){install.packages(new.packages)}

# On charge les packages nécessaires
library(data.table); library(terra)

# On charge la fonction qui permet de lire les sorties PHENOFIT#source(file.path(wd, "scripts", "functions", "read_phenofit_output.R"))
source("phenofit/imfusion/scripts/functions/read_phenofit_output.R")

## 1. Lire les sorties pour l'espèce
species <- "fagus_sylvatica_new"

#'abies_alba
#'fagus_sylvatica
#'quercus_ilex
#'quercus_robur
#'fagus_sylvatica_new



sim.path <- file.path("phenofit","imfusion", "data", "simulations", "historical", "ERA5-LAND", "phenofit", species)

# on utilise la fonction chargée pour lire la sortie (Fitness ici)
#fitness <- read_phenofit_output(sim.path, years = c(1970:2000), 
                                #output.var = "Fitness")
fitness <- read_phenofit_output(sim.path, years = c(1985:2020), 
                                  output.var = "Fitness")


# on peut la transformer en raster (pour visualiser rapidement)
fitness.r <- rast(fitness, crs = "EPSG:4326")
plot(fitness.r)
par(mar = c(5, 4, 4, 2) - 1.2)
title(paste(species,'fitness'))


#calcul best threshold -> il faut d'abord extraire les valeurs sur les points sinon j'ai pas les données obs et je peux pas faire tss
#et il faut que je les récup avec les valeurs direct pas 0/1 comme sur l'autre script
#donc let's remettre l'atre bout de sript au dessus pour avoir une selection de points et extraire direct les valeurs 

# libraries
library(tidyr)
library(dplyr)
#library(rstan)
library(ggplot2)
library(sf)
library(pROC)
library(raster)
library(sp)


#' Liste des espèces à étudier
#' Abies alba
#' Fagus sylvatica
#' (Picea abies , Pinus sylvestris)
#' Quercus ilex
#' Quercus robur

sp="Fagus sylvatica"

#choix du bon modele
liste_modele <- read.csv("model_fit_safetymargin/output_safmarg_era.csv", sep = ";", header = TRUE)

liste_mod_filt<- liste_modele |>
  filter(species==sp)

mod <- case_when(
  liste_mod_filt$rhat[1] < 1.5 ~ "_2sm",
  liste_mod_filt$rhat[2] < 1.5 ~ "_fsm",
  liste_mod_filt$rhat[3] < 1.5 ~ "_hsm",
  TRUE ~ "_none")
print(mod)

#appeler le modèle
load(paste0("model_fit_safetymargin/",sp,mod,".RData")) #appelle les données 
sp.trait=readRDS("data/df.species") |> 
  filter(species==sp)

#appeler le modèle fsm
#load(paste0("model_fit_safetymargin/",sp,"_fsm.RData")) #appelle les données 
#sp.trait=readRDS("data/df.species") |> 
# filter(species==sp)


occurence<-readRDS("data/occurence") |> 
  filter(species==sp) |> #filtre la bonne espèce dans "occurence'
  mutate(hsm=(psi_eraday_real/1000)-sp.trait$px,
         fsm=tmin_era-sp.trait$lt50) 
#' mutate remplace la valeur contenue dans hsm et fsm par leur nouvelle valeur calculée 
#' ici en fonction des données climatiques d'entrée (psi et tmin)


# get pars posteriors
post<-as.data.frame(fit.2var)

post |> 
  pivot_longer(cols=everything()) |> #inverse les lignes et les colonnes 
  ggplot(aes(value))+ #crée juste le plot de base / la première couche 
  geom_density()+
  #'crée le graphe de densité de probabilité des data mais on voit rien car
  #'les échelles de chaque variables sont trop différentes 
  facet_wrap(~name,ncol=1,scales="free") 
#'permet de faire un graphe pour chaque variable que l'on souhaite 
#'ici on appelle "name" qui doit être une liste avec toutes les variables K,r,t,
#'lp (? c'est quoi lp) pour hsm et fsm , donc on obtient 1 graphe par variable 

post |> 
  mutate(nsamp=row_number()) |> #attribue un num unique à chaque ligne -> ici 1 
  #'ligne = 1 combinaison des variables K,r,t ... 
  pivot_longer(cols=-nsamp) |> #passage en lignes sauf nsamp (on obtient à la suite K1,r1,t1,K2,r2,t2 ...)
  ggplot(aes(nsamp,value))+ #aes = aesthetic
  geom_line()+ #la contrairement à geomdensity , il prend juste les points un a un du "tirage" qui a été fait 
  facet_wrap(~name,ncol=1,scales="free_y")  


# plot map of Fagus presence 
worldmap <- st_as_sf(rworldmap::getMap(resolution = "high")) #installer rworldmap et rworldxtra pour faire tourner
#c'est litteralement la carte du monde : code des pays 

occurence |>
  group_by(presence) |> #regroupe  les données en f de la colonne presence
  sample_n(1000) #pioche 1 000 echantillons 
  #ggplot()+
  #geom_point(aes(x=x,y=y,color=as.factor(presence)))+ #les points en couleur
  #geom_sf(data=worldmap,fill=NA)+ #on affiche la carte dessous sans couleur
  #xlim(-20,35)+
  #ylim(27,73)  
  #'sort la carte du monde (worldmap), xlim et ylim pour 
  #'les coordonnées limites de la carte, coloré en f de "présence"
  
  # plot safety margin
  occurence |> 
  filter(!is.na(hsm),!is.na(fsm)) |> #filtre en enlevant les lignes n'ayant pas de hsm ou fsm 
  mutate(hsm=cut(hsm,
                 breaks=c(-Inf, -15, -10, -5,-3, -2, -1, 0,1,2,3,5,Inf)), #de - l'infini à +l'infini
         fsm=cut(fsm,
                 breaks=c(-Inf, -50,-40, -30,-20,-15,-10,-5,0,5,10,Inf))
  ) |>
  #'crée des classes d'intervalles avec cut (il faut couper) et breaks (points de coupure) de hsm et fsm
  sample_n(20000) |>
  ggplot(aes(x=x,y=y,color=hsm))+
  geom_point()+
  scale_color_brewer(palette="RdYlBu")
#sort une carte des hsm avec une couleur par classe (ici il y a tellement de 
#'points qu'on voit l'Europe sans même faire appel à la carte)

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
  scale_color_brewer(palette="RdYlBu") #même chose mais pour fsm 



# trajectory with uncertainty
fsm.05<- quantile(occurence$fsm,prob=0.05,na.rm=TRUE)[[1]] 
fsm.95<- quantile(occurence$fsm,prob=0.95,na.rm=TRUE)[[1]]
hsm.05<- quantile(occurence$hsm,prob=0.05,na.rm=TRUE)[[1]]
hsm.95<- quantile(occurence$hsm,prob=0.95,na.rm=TRUE)[[1]]


#pour tirer des valeurs en excluant les extrêmes, on prend entre 0,05 et 0,95

post |> 
  crossing(data.frame(xsm_val=c(seq(fsm.05,
                                    fsm.95,
                                    length.out=100), #on tire 100 fsm et 100 hsm
                                seq(hsm.05,
                                    hsm.95,#pou rne pas être limité par la hsm
                                    length.out=100)),# on les combine avec crossing :
                      #toutes les valeurs tirées avec tous les jeux de variables
                      xsm_name=c(rep("fsm",100),
                                 rep("hsm",100)),
                      xsm_name_l=c(rep("Frost safety margins (°C)",100),
                                   rep("Hydraulic safety margins (MPa)",100)))) |> 
  mutate(pred=case_when(xsm_name=="fsm"~K_int/((1+exp(-r_hsm*(hsm.95-t_hsm)))*    
                                                 (1+exp(-r_fsm*(xsm_val-t_fsm)))),
                        xsm_name=="hsm"~K_int/((1+exp(-r_fsm*(fsm.95-t_fsm)))*
                                                 (1+exp(-r_hsm*(xsm_val-t_hsm))))))|> 
  
  group_by(xsm_val,xsm_name_l) |> 
  summarise(med_pred=median(pred),
            q05_pred=quantile(pred,prob=0.05),
            q95_pred=quantile(pred,prob=0.95)) |> #pour avoir un intervalle de confiance
  
  ggplot()+
  geom_line(aes(xsm_val,med_pred))+
  geom_ribbon(aes(x=xsm_val,ymin=q05_pred,ymax=q95_pred),alpha=0.2)+ #fait apparitre la zone grisée/intervalle de c
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
post_med<-apply(post,MARGIN=2,median) #on prend la médiane de post
tss<-occurence |> 
  mutate(pred=post_med[["K_int"]]/
           ((1+exp(-post_med[["r_fsm"]]*(fsm-post_med[["t_fsm"]])))*
              (1+exp(-post_med[["r_hsm"]]*(hsm-post_med[["t_hsm"]]))))) |> 
  filter(!is.na(fsm),!is.na(hsm))
#on lui donne occurence et on calcule les proba de présence 

calc_tss <- function(threshold, observed, predicted_probs) {
  observed=as.factor(observed) #change les données observées en facteurs 
  predicted <- ifelse(predicted_probs > threshold, 1, 0)
  predicted=factor(predicted,levels=c(0,1)) #change en facteurs ordonnés 0 et 1 
  conf_matrix <- table(observed, predicted) 
  TP <- conf_matrix[2, 2] #associe la valeur dans le tableau à une variable TP -> predits+obs
  FN <- conf_matrix[2, 1] #pas pred mais obs
  TN <- conf_matrix[1, 1] #pas pred pas obs
  FP <- conf_matrix[1, 2] #pred mais pas obs
  sensitivity <- TP / (TP + FN) #les  perfectsobs sur tous les obs
  specificity <- TN / (TN + FP) #les perfectspasobs sur tous les pas obs
  TSS <- sensitivity + specificity - 1
  return(TSS)
}
observed <- tss$presence
pred<- tss$pred
threshold_max<-post_med[["K_int"]]
thresholds <- seq(0,threshold_max, length.out=100) #tirage de ceux sur lesquels on va appliquer la fonction
tss_sfm <- sapply(thresholds, calc_tss, observed, pred)


#on cherche a connaitre le threshold qui maximise le TSS

opt_prob=thresholds[which.max(tss_sfm)]


# Utilisation de sample pour sélectionner 5000 indices aléatoires parmi les indices des points existants
indices_aleatoires <- sample(length(tss$x),241045)#

# Sélection des points correspondant à ces indices
x_aleatoire <- tss$x[indices_aleatoires]
y_aleatoire <- tss$y[indices_aleatoires]

# Affichage des points générés

ggplot()+
  geom_point(aes(x=x_aleatoire,y=y_aleatoire))+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)

tss_ech <- tss |>
  filter(x%in%x_aleatoire)

tss_ech$pred <- ifelse(
  tss_ech$pred < opt_prob, # si la fitness est inférieure au seuil:
  0, # on considère que l'espèce est absente (=0)
  1 # sinon: on considère que l'espèce est présente (=1)
)

tss_ech |> 
  #mutate(pred_presence=pred>opt_prob) |> en gros si je fais tout apparaitre sur la carte pres
  #et absence les abs "cachent" les pres, détail d'affichage dc ok pour là 
  group_by(presence) |> 
  #sample_n(600) |> 
  pivot_longer(cols=c("presence","pred")) |> 
  ggplot()+
  geom_point(aes(x=x,y=y,color=value))+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)+
  facet_wrap(~name)+
  ggtitle(paste0("Prediction vs observation de ",sp))

#ok donc là je visualise les points, il faut maintenant extraires ces valeurs du raster

# extraction du raster sur les points cibles

points_alea_df <- data.frame(x = tss_ech$x, y = tss_ech$y)

pred_phenofit_df <- extract(fitness.r,points_alea_df)  #bien avoir chargé library sp et raster
tss_ech<- mutate(tss_ech,pred_phenofit=pred_phenofit_df$val)


#calcul du treshold  

calc_tss <- function(threshold, observed, predicted_probs) {
  observed=as.factor(observed) #change les données observées en facteurs 
  predicted <- ifelse(predicted_probs > threshold, 1, 0)
  predicted=factor(predicted,levels=c(0,1)) #change en facteurs ordonnés 0 et 1 
  conf_matrix <- table(observed, predicted) 
  TP <- conf_matrix[2, 2] #associe la valeur dans le tableau à une variable TP -> predits+obs
  FN <- conf_matrix[2, 1] #pas pred mais obs
  TN <- conf_matrix[1, 1] #pas pred pas obs
  FP <- conf_matrix[1, 2] #pred mais pas obs
  sensitivity <- TP / (TP + FN) #les  perfectsobs sur tous les obs
  specificity <- TN / (TN + FP) #les perfectspasobs sur tous les pas obs
  TSS <- sensitivity + specificity - 1
  return(TSS)
}
observed <- tss_ech$presence
pred_ph<- tss_ech$pred_phenofit
threshold_max_ph<-max(tss_ech$pred_phenofit) #on prend le max pour phenofit (ou 1 c'est pareil de toute façon)
thresholds_ph <- seq(0,threshold_max_ph, length.out=100) 
tss_sfm_ph <- sapply(thresholds_ph, calc_tss, observed, pred_ph)
#on cherche a connaitre le threshold qui maximise le TSS

opt_prob_ph=thresholds_ph[which.max(tss_sfm_ph)]
print(opt_prob_ph)



#remplacer par 0 et 1 et sortir les cartes 

tss_ech$pred_phenofit <- ifelse(
  tss_ech$pred_phenofit < opt_prob_ph, # si la fitness est inférieure au seuil:
  0, # on considère que l'espèce est absente (=0)
  1 # sinon: on considère que l'espèce est présente (=1)
)


tss_ech |> 
  #mutate(pred_presence=pred>opt_prob) |> #presence vraie ou fausse
  #mutate(pred_presence_ph=pred_phenofit>opt_prob_ph)
  group_by(presence) |> 
  #sample_n(600) |> #il veut 693 car c'est le nb de 1 présence 
  pivot_longer(cols=c("presence","pred","pred_phenofit")) |> 
  #subset(value == 1) |> #pour n'afficher que les présences 
  ggplot()+
  geom_point(aes(x=x,y=y,color=value,alpha= ifelse(value == 0, 0.5, 1)))+
  scale_color_gradient(low = "lightgreen", high = "darkgreen")+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)+
  facet_wrap(~name)+
  ggtitle(paste0("Prediction vs observation MA et pheno de ",sp)) 


##calcul de la similarité avec la présence réelle 

tss_ech$simil_ph<- ifelse(tss_ech$pred_phenofit==tss$presence,1,0)

similitude_ph <- sum(tss_ech$simil==1)/nrow(tss_ech) * 100

print(similitude_ph) 

tss_ech$simil<- ifelse(tss_ech$pred==tss$presence,1,0)

similitude <- sum(tss_ech$simil==1)/nrow(tss_ech) * 100

print(similitude) 

#'similitude plus grande avec MA alors qu'on dirait pas forcément sur la carte
#'surement bcp plsu de 1 dans pheno alos qu'en présence réelle bcp d'absence aussi 
#'



  













