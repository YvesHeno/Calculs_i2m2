main <- function(Fichier_entree){
#import des fonctions de  calc consult
source("R/I2M2_v1.0.6_calc_consultYM.r")
print(paste("fichier entree :",Fichier_entree))
#version/indic
indic  <- "I2M2"
vIndic <- "v1.0.6"

## CHARGEMENT DES PACKAGES ----
dependencies <- c("dplyr", "sfsmisc", "tidyr")

loadDependencies <- function(dependencies) {
  suppressAll <- function(expr) {
    suppressPackageStartupMessages(suppressWarnings(expr))
  }

  lapply(dependencies,
         function(x)
         {
           suppressAll(library(x, character.only = TRUE))
         }
  )
  invisible()
}

loadDependencies(dependencies)

## IMPORT DES FICHIERS DE CONFIGURATION ----
Transcodage <- read.csv2("data/I2M2_params_transcodage.csv",
                         header = TRUE, stringsAsFactors = FALSE,
                         colClasses = c(CODE_TAXON   = "character",
                                        CODE_METHODE = "character"))
Base      <- read.csv2("data/I2M2_params_base.csv",
                       colClasses = c(cd_taxon = "character"))
Typo      <- read.csv2("data/I2M2_params_typo_I2M2.csv",
                       stringsAsFactors = FALSE)
Best      <- read.csv2("data/I2M2_params_best.csv")
Worst     <- read.csv2("data/I2M2_params_worst.csv")
De        <- read.csv2("data/I2M2_params_de.csv",
                       stringsAsFactors = FALSE)






## INITIALISATION DU TRAITEMENT ----
# Ne pas afficher les messages d'avis ni d'erreur
options(warn = -1)

# Fichier_entree           <- choose.files(caption = "Choisir le fichier", multi = FALSE,
#                                filters = cbind("Fichier texte (*.txt)", "*.txt"))
complementaire <- TRUE

# Initialisation de l'heure
heure_debut <- Sys.time()

##  IMPORT DES FICHIERS ----
# Import du fichier d'entree
data_entree <- read.table(Fichier_entree, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "\"",
                          colClasses = c(CODE_OPERATION = "character",
                                         CODE_STATION   = "character",
                                         CODE_TAXON     = "character")) %>%
  filter(RESULTAT > 0)

## INITIALISATION DU FICHIER DE SORTIE ----
paramsOut <- data.frame(CODE_PAR = c(8058, 8057, 8056, 8055, 8054, 7613, 8050),
                        LIB_PAR  = c("IndiceShannonI2M2",
                                     "AverageScorePerTaxonI2M2",
                                     "PolyvoltinismeI2M2",
                                     "OvovivipariteI2M2",
                                     "RichesseI2M2",
                                     "Ind Invert Multimetrique",
                                     "NbTaxonsI2M2Contributifs"),
                        stringsAsFactors = FALSE)

data_sortie <- funSortie(data_entree = data_entree,
                         paramsOut   = paramsOut,
                         CODE_OPERATION, CODE_STATION, DATE) %>%
  mutate(CODE_OPERATION = as.character(CODE_OPERATION),
         CODE_STATION   = as.character(CODE_STATION),
         DATE           = as.character(DATE))

## TAXONS CONTRIBUTIFS ----
TaxaContributifs <- funContributif(Table       = data_entree,
                                   Transcodage = Transcodage)

## CALCUL DE L'INDICE ----
metriques <- data_entree                                 %>%
  funHarmonisation(Table = ., Transcodage = Transcodage) %>%
  funEffBocaux(Table = .)                                %>%
  funMetriques(Table = ., Base = Base)

EQR <- metriques                  %>%
  funTypo(Table = ., Typo = Typo) %>%
  funEQR(Table = ., Worst = Worst, Best = Best)

I2M2 <- funI2M2(Table = EQR, De = De)

## RESULTATS COMPLEMENTAIRES ----
if (complementaire) {
  data_complementaire <- mutate(metriques,
                                SHANNON     = funArrondi(SHANNON, 4),
                                ASPT        = funArrondi(ASPT, 4),
                                POLYVOLTIN  = funArrondi(POLYVOLTIN, 4),
                                OVOVIVIPARE = funArrondi(OVOVIVIPARE, 4))
} else {
  data_complementaire <- NULL
}

## SORTIE DES RESULTATS ----
resultats <- bind_rows(
  funFormat(EQR),
  funFormat(I2M2),
  funFormat(TaxaContributifs)
)

data_sortie <- left_join(x  = data_sortie,
                         y  = resultats,
                         by = c("CODE_OPERATION", "CODE_PAR"))


Commentaires <- funCommentaires(Table            = data_sortie,
                                TaxaContributifs = TaxaContributifs)

data_sortie <- left_join(x = data_sortie,
                         y = Commentaires,
                         by = c("CODE_OPERATION", "CODE_PAR")) %>%
  mutate(RESULTAT = ifelse(CODE_PAR == 8050,
                           RESULTAT,
                           funArrondi(RESULTAT, 4)),
         COMMENTAIRE = ifelse(is.na(COMMENTAIRE), "", COMMENTAIRE))

#le fichier de sortie a pour debut de nom le fichier d'netree, on enleve le .txt
#pour ne pas avoir de .txt au milieu
# de plus la liste des fichiers d'entree regarde tous les fichiers qui ont
#.txt dans leur nom
Ficentreemod <- str_remove(Fichier_entree,".txt")
#fichierResultat               <- paste0("data/",indic, "_", vIndic, "_resultats.csv")
fichierResultat               <- paste0(Ficentreemod, "-",vIndic, "_resultats.csv")
#fichierResultatComplementaire <- paste0("data",indic, "_", vIndic,"_resultats_complementaires.csv")
fichierResultatComplementaire <- paste0(Ficentreemod, "_", vIndic,"_resultats_complementaires.csv")
funResult(indic               = indic,
          vIndic              = vIndic,
          heure_debut         = heure_debut,
          data_sortie         = data_sortie,
          data_complementaire = data_complementaire,
          complementaire      = complementaire,
          file                = fichierResultat,
          file_complementaire = fichierResultatComplementaire)

return(data_sortie)
print(paste("fin du traitement, résultats dans ",fichierResultat))
print("**********************************************************")
}
sortie <- main("data/I2M2_entree_op100.txt")

Liste_fichiers <-  list.files("data/",pattern="*.txt",full.names=T)
test <- map_df(.x=Liste_fichiers,
    .f=main)

#test2 <- test %>% reduce(rbind)
Extracti2m2<-filter(test,CODE_PAR=='7613') %>% distinct() #distinct pour enlever les doublons
colnames(Extracti2m2)[2] <- "cdstation"

stations <- "https://geobretagne.fr/geoserver/dreal_b/stations_qualite_eau_surf/wfs?SERVICE=WFS&REQUEST=GetCapabilities"
leaflet::leaflet() %>% setView(lng = -3, lat = 48, zoom = 8) %>%
  leaflet::addWMSTiles(stations, layers = "stations_qualite_eau_surf")


st <- leaflet::leaflet() %>% setView(lng = -3, lat = 48, zoom = 8) %>%
  addWMSTiles(stations, layers = "stations_qualite_eau_surf")

stations <- sf::st_read("couches/stations_qualite_eau_surf.shp",options = "ENCODING=WINDOWS-1252") %>% filter(y>687000)# il y a une station en
#afrique, on la vire
library(sf)
#mapview::mapview(stations)
indicesgeom <-left_join(Extracti2m2,stations,by="cdstation") %>% st_as_sf()
mapview::mapview(indicesgeom)

