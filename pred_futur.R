# libraries
library(tidyr)
library(dplyr)
#library(rstan)
library(ggplot2)
library(sf)
library(pROC)
library(raster)
library(sp)
library(data.table)
library(terra)

library(ncdf4)


#ouvrir le fichier des températures 

carte_temp <- "clim_futur/frost_future/IPSL/ssp245/tasminAdjust_day_IPSL-CM6A-LR_ssp245_r1i1p1f1_gr010_TCDF-CDFT23-ERA5Land-1981-2010_2046-2055_percentile2.nc"

#2036-2045
#2046-2055
#2056-2065
#2066-2075
#2076-2085
#2086-2100


carte_temp_rast <- raster(carte_temp)



# Afficher les valeurs min et max de vos coordonnées géographiques avant conversion
print(extent(carte_temp_rast))

# Ajuster les coordonnées géographiques
extent(carte_temp_rast) <- extent(-180, 180, -90, 90)

# Afficher les valeurs min et max de vos coordonnées géographiques après conversion
print(extent(carte_temp_rast))


# Définir l'étendue de l'Europe est
europe_extent_east <- extent(-180, -140, 35, 70)  # xmin, xmax, ymin, ymax

# Ouvrir une nouvelle fenêtre pour la figure
windows()

# Visualiser le raster avec la nouvelle étendue
plot(carte_temp_rast, ext=europe_extent_east)


# Créer un nouveau raster vide avec la même résolution et la même étendue que le raster d'origine
europe_raster_east <- raster(extent(europe_extent_east), res = res(carte_temp_rast))

# Remplir le nouveau raster avec les valeurs du raster d'origine uniquement dans la partie est de l'Europe
europe_raster_east[] <- crop(carte_temp_rast, europe_extent_east)[]

# Afficher le nouveau raster
windows()
plot(europe_raster_east)



# Définir l'étendue de l'Europe à l'ouest
europe_extent_west <- extent(160, 180, 35, 70)  # xmin, xmax, ymin, ymax


# Visualiser le raster avec la nouvelle étendue
windows()
plot(carte_temp_rast, ext=europe_extent_west)


# Créer un nouveau raster vide avec la même résolution et la même étendue que le raster d'origine
europe_raster_west <- raster(extent(europe_extent_west), res = res(carte_temp_rast))

# Remplir le nouveau raster avec les valeurs du raster d'origine uniquement dans la partie est de l'Europe
europe_raster_west[] <- crop(carte_temp_rast, europe_extent_west)[]

# Afficher le nouveau raster
windows()
plot(europe_raster_west)

------------------------------------------------------------------------------------------------------
  
#Redimentionner les nouveaux raster 

extent(europe_raster_east) <- extent(0, 40, 35, 70)
extent(europe_raster_west) <- extent(-20,0,35,70)

# Afficher la carte 

europe_raster_combined <- merge(europe_raster_east, europe_raster_west)
windows()
plot(europe_raster_combined)


--------------------------------------------------------------------------------------------------

####Utiliser les deux raster pour extraire les valeurs et les mettre dans une même colonne 
##modèle anne 


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
  sample_n(1000) |> #pioche 1 000 echantillons 
  ggplot()+
  geom_point(aes(x=x,y=y,color=as.factor(presence)))+ #les points en couleur
  geom_sf(data=worldmap,fill=NA)+ #on affiche la carte dessous sans couleur
  xlim(-20,35)+
  ylim(27,73)  
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
#'sort une carte des hsm avec une couleur par classe (ici il y a tellement de 
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
hsm.999<-quantile(occurence$hsm,prob=0.999,na.rm=TRUE)[[1]] 

#pour tirer des valeurs en excluant les extrêmes, on prend entre 0,05 et 0,95

post |> 
  crossing(data.frame(xsm_val=c(seq(fsm.05,
                                    fsm.95,
                                    length.out=100), #on tire 100 fsm et 100 hsm
                                seq(hsm.05,
                                    hsm.95,
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
              (1+exp(-post_med[["r_hsm"]]*(hsm.999-post_med[["t_hsm"]]))))) |> #pour fixer hsm
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



## make prediction and plot

tss |> 
  mutate(pred_presence=pred>opt_prob) |> 
  group_by(presence) |> 
  #sample_n(600) |> 
  pivot_longer(cols=c("presence","pred_presence")) |> 
  ggplot()+
  geom_point(aes(x=x,y=y,color=value))+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)+
  facet_wrap(~name)+
  ggtitle(paste0("Prediction vs observation de ",sp, " modele ",mod," hsm fixe"))

--------------------------------------------------------------------------------------------------------------

##extraire les valeurs des rasters 

# Utilisation de sample pour sélectionner 5000 indices aléatoires parmi les indices des points existants
indices_aleatoires <- sample(length(occurence$x),241045) #241045

# Sélection des points correspondant à ces indices
x_aleatoire <- occurence$x[indices_aleatoires]
y_aleatoire <- occurence$y[indices_aleatoires]

# Affichage des points générés

ggplot()+
  geom_point(aes(x=x_aleatoire,y=y_aleatoire))+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)

occurence_ech_futur <- occurence |>
  filter(x%in%x_aleatoire)


# extraction du raster sur les points cibles

points_extract <- data.frame(x = occurence_ech_futur$x, y = occurence_ech_futur$y)

temp_futur_list <- extract(europe_raster_combined,points_extract)  #bien avoir chargé library sp et raster

# Convertir la liste de valeurs en un vecteur
temp_futur <- unlist(temp_futur_list)

# Ajouter les valeurs extraites à votre dataframe initial
occurence_ech_futur$tmin_futur <- temp_futur

occurence_ech_futur$tmin_futur <- occurence_ech_futur$tmin_futur - 273.15


-----------------------------------------------------------------------------------------------------

#modèle avec nouvelles temp

occurence_ech_futur |> 
  mutate(hsm=(psi_eraday_real/1000)-sp.trait$px,
         fsm=tmin_futur-sp.trait$lt50) 
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
  sample_n(1000) |> #pioche 1 000 echantillons 
  ggplot()+
  geom_point(aes(x=x,y=y,color=as.factor(presence)))+ #les points en couleur
  geom_sf(data=worldmap,fill=NA)+ #on affiche la carte dessous sans couleur
  xlim(-20,35)+
  ylim(27,73)  
#'sort la carte du monde (worldmap), xlim et ylim pour 
#'les coordonnées limites de la carte, coloré en f de "présence"

# plot safety margin
occurence_ech_futur |> 
  filter(!is.na(hsm),!is.na(fsm)) |> #filtre en enlevant les lignes n'ayant pas de hsm ou fsm 
  mutate(hsm=cut(hsm,
                 breaks=c(-Inf, -15, -10, -5,-3, -2, -1, 0,1,2,3,5,Inf)), #de - l'infini à +l'infini
         fsm=cut(fsm,
                 breaks=c(-Inf, -50,-40, -30,-20,-15,-10,-5,0,5,10,Inf))
  ) |> 
  #'crée des classes d'intervalles avec cut (il faut couper) et breaks (points de coupure) de hsm et fsm
  sample_n(20000) |> #-------------------------------------------------------------------------------------------ligne changée 20000-> 2000 car ech pour faire les tests + vite
  ggplot(aes(x=x,y=y,color=hsm))+
  geom_point()+
  scale_color_brewer(palette="RdYlBu")
#'sort une carte des hsm avec une couleur par classe (ici il y a tellement de 
#'points qu'on voit l'Europe sans même faire appel à la carte)

occurence_ech_futur |> 
  mutate(hsm=(psi_eraday_real/1000)-sp.trait$px,
         fsm=tmin_futur-sp.trait$lt50) |> 
  filter(!is.na(hsm),!is.na(fsm)) |> 
  mutate(hsm=cut(hsm,
                 breaks=c(-Inf, -15, -10, -5,-3, -2, -1, 0,1,2,3,5,Inf)),
         fsm=cut(fsm,
                 breaks=c(-Inf, -50,-40, -30,-20,-15,-10,-5,0,5,10,Inf))
  ) |> 
  sample_n(20000) |> #--------------------------------------------------------------même chose
  ggplot(aes(x=x,y=y,color=fsm))+
  geom_point()+
  scale_color_brewer(palette="RdYlBu") #même chose mais pour fsm 



# trajectory with uncertainty
fsm.05<- quantile(occurence$fsm,prob=0.05,na.rm=TRUE)[[1]] 
fsm.95<- quantile(occurence$fsm,prob=0.95,na.rm=TRUE)[[1]]
hsm.05<- quantile(occurence$hsm,prob=0.05,na.rm=TRUE)[[1]]
hsm.95<- quantile(occurence$hsm,prob=0.95,na.rm=TRUE)[[1]]
hsm.999<-quantile(occurence$hsm,prob=0.999,na.rm=TRUE)[[1]] 

#pour tirer des valeurs en excluant les extrêmes, on prend entre 0,05 et 0,95

post |> 
  crossing(data.frame(xsm_val=c(seq(fsm.05,
                                    fsm.95,
                                    length.out=100), #on tire 100 fsm et 100 hsm
                                seq(hsm.05,
                                    hsm.95,
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
tss_ech_futur<-occurence_ech_futur |> 
  mutate(pred=post_med[["K_int"]]/
           ((1+exp(-post_med[["r_fsm"]]*(fsm-post_med[["t_fsm"]])))*
              (1+exp(-post_med[["r_hsm"]]*(hsm.999-post_med[["t_hsm"]]))))) |> #pour fixer hsm
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
observed <- tss_ech_futur$presence
pred<- tss_ech_futur$pred
threshold_max<-post_med[["K_int"]]
thresholds <- seq(0,threshold_max, length.out=100) #tirage de ceux sur lesquels on va appliquer la fonction
tss_sfm <- sapply(thresholds, calc_tss, observed, pred)
#on cherche a connaitre le threshold qui maximise le TSS

opt_prob=thresholds[which.max(tss_sfm)]



## make prediction and plot

tss_ech_futur |> 
  mutate(pred_presence=pred>opt_prob) |> 
  group_by(presence) |> 
  pivot_longer(cols=c("presence","pred_presence")) |> 
  ggplot()+
  geom_point(aes(x=x,y=y,color=value,alpha= ifelse(value == 0, 0.5, 1)))+
  scale_color_gradient(low = "lightgreen", high = "darkgreen")+
  geom_sf(data=worldmap,fill=NA)+
  xlim(-20,35)+
  ylim(27,73)+
  facet_wrap(~name)+
  ggtitle(paste0("Prediction vs observation de ",sp, " modele ",mod," hsm fixe", "futur")) #automatiser ensuite pour toutes les plages de dates


## Dans l'idée j'aimerai faire un facet wrap avec toutes les cartes  qui sortent en même temps. Et je pense enlever la présence actuelle car ce serait pas pertinent 





