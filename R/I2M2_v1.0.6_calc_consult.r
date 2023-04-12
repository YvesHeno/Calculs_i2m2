# Type d'algorithme : I2M2
# Auteur(s)         : Cedric MONDY, Delphine CORNEIL
# Date              : 2019-11-12
# Version           : 1.0.6
# Interpreteur	   	: R version 3.5.3 (2019-03-11)
# Pre-requis        : Packages dplyr, sfsmisc, tidyr
# Fichiers lies   	: I2M2_params_base.csv, I2M2_params_best.csv, I2M2_params_de.csv, I2M2_params_transcodage.csv, I2M2_params_typo_I2M2.csv, I2M2_params_worst.csv
# Commentaires 	  	: Mondy CP, Villeneuve B, Archaimbault V, Usseglio-Polatera
# P. (2012) A new macroinvertebrate-based multimetric index (I2M2) to evaluate
# ecological quality of French wadeable streams fulfilling the WFD demands: A
# taxonomical and trait approach. Ecological indicators, 18: 452-67.
# Usseglio-Polatera, P. & Mondy, C. (2011) Developpement et optimisation de
# l'indice biologique macroinvertebres benthiques (I2M2) pour les cours d'eau.
# Partenariat Onema / UPV-Metz - LIEBE - UMR-CNRS 7146, 27p.

# Copyright 2018 Cedric MONDY, Delphine CORNEIL
# Ce programme est un logiciel libre; vous pouvez le redistribuer ou le modifier
# suivant les termes de la GNU General Public License telle que publiee par la
# Free Software Foundation; soit la version 3 de la licence, soit (a votre gre)
# toute version ulterieure.
# Ce programme est distribue dans l'espoir qu'il sera utile, mais SANS AUCUNE
# GARANTIE; sans meme la garantie tacite de QUALITE MARCHANDE ou d'ADEQUATION A
# UN BUT PARTICULIER. Consultez la GNU General Public License pour plus de
# details.
# Vous devez avoir recu une copie de la GNU General Public License en meme temps
# que ce programme; si ce n'est pas le cas, consultez
# <http://www.gnu.org/licenses>.

## VERSION ----
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
Transcodage <- read.csv2("I2M2_params_transcodage.csv", 
                         header = TRUE, stringsAsFactors = FALSE, 
                         colClasses = c(CODE_TAXON   = "character",
                                        CODE_METHODE = "character"))
Base      <- read.csv2("I2M2_params_base.csv",
                       colClasses = c(cd_taxon = "character"))
Typo      <- read.csv2("I2M2_params_typo_I2M2.csv",
                       stringsAsFactors = FALSE)
Best      <- read.csv2("I2M2_params_best.csv")
Worst     <- read.csv2("I2M2_params_worst.csv")
De        <- read.csv2("I2M2_params_de.csv",
                       stringsAsFactors = FALSE)

## DECLARATION DES FONCTIONS ----
## Fonction determinant le nombre de taxons contributifs ou non
funContributif <- function(Table, Transcodage) {
  left_join(x = Table, y = Transcodage, by = "CODE_TAXON") %>%
    group_by(CODE_OPERATION)                               %>%
    summarise(nbTotal          = n_distinct(CODE_TAXON),
              nbContributif    = n_distinct(CODE_TAXON[! is.na(CODE_METHODE)]),
              nbNonContributif = n_distinct(CODE_TAXON[is.na(CODE_METHODE)]),
              taxaNonContributifs = CODE_TAXON[is.na(CODE_METHODE)] %>%
                unique()                                            %>%
                paste(collapse = ", "))
  
}

## Fonction permettant d'harmoniser les listes taxonomiques
funHarmonisation <- function(Table, Transcodage) {
  
  getAllTaxa <- function(df) {
    summ <- group_by(df, CODE_OPERATION) %>%
      summarise(taxa = paste(unique(CODE_TAXON), collapse = ", "))
    
    left_join(df, summ, by = "CODE_OPERATION")      %>%
      split(x = ., f = .$CODE_OPERATION)            %>%
      lapply(function(df.i) {
        expand.grid(CODE_OPERATION = unique(df.i$CODE_OPERATION),
                    CODE_PHASE     = unique(df.i$CODE_PHASE), 
                    TYPO_NATIONALE = unique(df.i$TYPO_NATIONALE),
                    CODE_TAXON     = strsplit(df.i$taxa, split = ", ")[[1]],
                    stringsAsFactors = FALSE)
      })                                            %>%
      bind_rows()                                   %>%
      left_join(x = ., y = df, 
                by = c("CODE_OPERATION", "TYPO_NATIONALE",
                       "CODE_PHASE", "CODE_TAXON")) %>%
      mutate(RESULTAT = ifelse(is.na(RESULTAT), 0,
                               RESULTAT))           %>%
      as.tbl()
  }
  
  funRedistribution <- function(df) {
    if (nrow(filter(df, !is.na(REDISTRIBUTION))) == 0 |
        is.na(unique(df$FAMILLE))                     |
        sum(df$RESULTAT) == 0) {
      0
    } else {
      filter(df, !is.na(REDISTRIBUTION)) %>%
        '$'("RESULTAT")                %>%
        (function(x) {
          roundfixS(x * df$PROP)
        })
    }
  }
  
  Table_transformee <- mutate(Table, 
                              # Taxons en presence/absence
                              RESULTAT = ifelse(CODE_REMARQUE == 4, 
                                                1, RESULTAT))               %>%
    group_by(CODE_OPERATION, TYPO_NATIONALE,
             CODE_PHASE, CODE_TAXON)                   %>%
    summarise(RESULTAT       = sum(RESULTAT))          %>%
    getAllTaxa(df = .)                                 %>%
    left_join(Transcodage, by = "CODE_TAXON")          %>%
    filter(!is.na(CODE_METHODE))                       %>%
    # Regroupe les taxons determines a un niveau plus precis que demande
    group_by(CODE_OPERATION, TYPO_NATIONALE,
             CODE_PHASE, CODE_METHODE)                 %>%
    summarise(RESULTAT = sum(RESULTAT))                %>%
    left_join(select(Transcodage, -CODE_METHODE),
              by = c("CODE_METHODE" = "CODE_TAXON"))   %>%
    # Taxons en presence/absence
    mutate(RESULTAT = ifelse(QUANTIFICATION == "presence" &
                               RESULTAT > 0,
                             1, RESULTAT))             %>%
    ungroup()                                          %>%
    mutate(op_famille = paste0(CODE_OPERATION, "_", FAMILLE))
  
  # Redisribue les taxons determines a un niveau insuffisant
  ## liste les op x familles pour lesquelles doit redistribuer
  redistributions <- Table_transformee                       %>%
    group_by(op_famille)                                     %>%
    summarise(nb_redistri = length(na.omit(REDISTRIBUTION)),
              nb_non_redistri = sum(is.na(REDISTRIBUTION)))  %>%
    filter(nb_redistri > 0 & nb_non_redistri > 0)            %>%
    '$'("op_famille")
  
  Table_sans_redistribution <- filter(Table_transformee,
                                      ! op_famille %in% redistributions)
  Table_avec_redistribution <- filter(Table_transformee,
                                      op_famille %in% redistributions)
  
  Table_redistribuee <- Table_avec_redistribution %>%
    ## Calcul des proportions des taxons enfants
    left_join(x = ., 
              y = group_by(., CODE_OPERATION, FAMILLE) %>%
                filter(is.na(REDISTRIBUTION))          %>%
                summarise(TOTAL = sum(RESULTAT)),
              by = c("CODE_OPERATION", "FAMILLE")) %>%
    left_join(x = .,
              y = group_by(., CODE_OPERATION, FAMILLE,
                           CODE_METHODE) %>%
                summarise(PROP = sum(RESULTAT) / unique(TOTAL)),
              by = c("CODE_OPERATION", "CODE_METHODE",
                     "FAMILLE"))                   %>%
    mutate(PROP = ifelse(is.na(TOTAL), 0, 
                         ifelse(!is.na(REDISTRIBUTION), -1,
                                PROP)))            %>%
    ungroup(.)                                     %>%
    ## Redistribue les effectifs
    mutate(fac = paste0(CODE_OPERATION, "_", 
                        CODE_PHASE, "_", FAMILLE)) %>%
    (function(df) {
      split(df, df$fac) %>%
        lapply(FUN = function(df) {
          mutate(df, RESULTAT = RESULTAT + funRedistribution(df))
        })              %>%
        bind_rows()
    })
  
  if (nrow(Table_redistribuee) > 0) {
    bind_rows(
      select(Table_sans_redistribution,
             CODE_OPERATION, TYPO_NATIONALE, 
             CODE_PHASE, CODE_METHODE, RESULTAT),
      select(Table_redistribuee,
             CODE_OPERATION, TYPO_NATIONALE,
             CODE_PHASE, CODE_METHODE, RESULTAT)
    )                                                 %>%
      arrange(CODE_OPERATION, CODE_PHASE, CODE_METHODE) %>%
      filter(RESULTAT > 0)
    
  } else {
    select(Table_sans_redistribution,
           CODE_OPERATION, TYPO_NATIONALE, 
           CODE_PHASE, CODE_METHODE, RESULTAT) %>%
      arrange(CODE_OPERATION, CODE_PHASE, as.numeric(CODE_METHODE)) %>%
      filter(RESULTAT > 0)
  }
}

## Fonction permettant de calculer les effectifs par regroupement de bocaux
funEffBocaux <- function(Table) {
  ungroup(Table) %>%
    (function(df) {
      bind_rows(
        filter(df, CODE_PHASE %in% c("A", "B", 1, 2))         %>%
          mutate(CODE_PHASE = "AB")                           %>%
          group_by(CODE_OPERATION, TYPO_NATIONALE,
                   CODE_PHASE, CODE_METHODE)  %>%
          summarise(RESULTAT   = sum(RESULTAT)),
        filter(df, CODE_PHASE %in% c("B", "C", 2, 3))         %>%
          mutate(CODE_PHASE = "BC")                           %>%
          group_by(CODE_OPERATION, TYPO_NATIONALE,
                   CODE_PHASE, CODE_METHODE)  %>%
          summarise(RESULTAT   = sum(RESULTAT)),
        filter(df, CODE_PHASE %in% c("A", "B", "C", 1, 2, 3)) %>%
          mutate(CODE_PHASE = "ABC")                          %>%
          group_by(CODE_OPERATION, TYPO_NATIONALE,
                   CODE_PHASE, CODE_METHODE)  %>%
          summarise(RESULTAT   = sum(RESULTAT))
      )
    })           %>%
    filter(RESULTAT > 0)
}

## Fonction permettant de calculer les metriques
funMetriques <- function(Table, Base) {
  ## Fonction permettant le calcul de la diversite de Shannon
  funShannon <- function(Table) {
    filter(Table, 
           CODE_PHASE == "AB",
           RESULTAT > 0)                            %>%
      group_by(CODE_OPERATION, TYPO_NATIONALE)      %>%
      mutate(PROPORTION = RESULTAT / sum(RESULTAT)) %>%
      summarise(SHANNON = -sum(PROPORTION * log2(PROPORTION)))
  }
  
  ## Fonction permettant de calculer la metrique aspt
  funASPT <- function(Table, Base) {
    filter(Table, 
           CODE_PHASE == "BC",
           RESULTAT > 0)                                   %>%
      left_join(Base, by = c("CODE_METHODE" = "cd_taxon")) %>%
      group_by(CODE_OPERATION, TYPO_NATIONALE, BMWP.taxo)  %>%
      summarise(RESULTAT = sum(RESULTAT),
                SCORE    = unique(BMWP.Original.Score))    %>%
      filter(! is.na(SCORE))                               %>%
      group_by(CODE_OPERATION, TYPO_NATIONALE)             %>%
      summarise(BMWP = sum(SCORE),
                ASPT = mean(SCORE))                        %>%
      select(CODE_OPERATION, TYPO_NATIONALE, ASPT)
  }
  
  ## Fonction permettant de calculer la metrique polyvoltin
  funPolyvoltin <- function(Table, Base) {
    filter(Table,
           CODE_PHASE == "ABC",
           RESULTAT > 0)                                   %>%
      left_join(Base, by = c("CODE_METHODE" = "cd_taxon")) %>%
      group_by(CODE_OPERATION, TYPO_NATIONALE)             %>%
      summarise(POLYVOLTIN = sum(log1p(RESULTAT) * polyvoltine, 
                                 na.rm = TRUE) /
                  sum(log1p(RESULTAT), na.rm = TRUE))
    
  }
  
  ## Fonction permettant de calculer la metrique ovovivipare
  funOvovivipare <- function(Table, Base) {
    filter(Table, 
           CODE_PHASE == "ABC",
           RESULTAT > 0)                                   %>%
      left_join(Base, by = c("CODE_METHODE" = "cd_taxon")) %>%
      group_by(CODE_OPERATION, TYPO_NATIONALE)             %>%
      summarise(OVOVIVIPARE = sum(log1p(RESULTAT) * ovoviviparity,
                                  na.rm = TRUE) /
                  sum(log1p(RESULTAT), na.rm = TRUE))
  }
  
  ## Fonction permettant le calcul de la richesse taxonomique
  funRichesse <- function(Table) {
    filter(Table, 
           CODE_PHASE == "ABC",
           RESULTAT > 0)                       %>%
      filter(RESULTAT > 0)                     %>%
      group_by(CODE_OPERATION, TYPO_NATIONALE) %>%
      summarise(RICHESSE = n_distinct(CODE_METHODE))
  }
  
  ## Calcul des metriques
  funShannon(Table)                                       %>%
    full_join(funASPT(Table, Base),        
              by = c("CODE_OPERATION", "TYPO_NATIONALE")) %>%
    full_join(funPolyvoltin(Table, Base),  
              by = c("CODE_OPERATION", "TYPO_NATIONALE")) %>%
    full_join(funOvovivipare(Table, Base), 
              by = c("CODE_OPERATION", "TYPO_NATIONALE")) %>%
    full_join(funRichesse(Table),          
              by = c("CODE_OPERATION", "TYPO_NATIONALE"))
}

## Fonction permettant d'attribuer les typologies simplifiees
funTypo <- function(Table, Typo) {
  left_join(Table, 
            select(Typo, TYPO_NATIONALE, TYPO_I2M2),
            by = "TYPO_NATIONALE")
}

## Fonction permettant de calculer les EQR des metriques
funEQR <- function(Table, Worst, Best) {
  Table                                               %>%
    select(-TYPO_NATIONALE)                           %>%
    gather(key = METRIQUE, value = RESULTAT,
           -CODE_OPERATION, -TYPO_I2M2)               %>%
    left_join(gather(Worst, key = METRIQUE, value = WORST), 
              by = "METRIQUE")                        %>%
    left_join(gather(Best, key = METRIQUE, value = BEST, 
                     -TYPO_I2M2),
              by = c("TYPO_I2M2", "METRIQUE"))        %>%
    mutate(EQR = (RESULTAT - WORST) / (BEST - WORST)) %>%
    mutate(EQR = ifelse(EQR < 0, 0, 
                        ifelse(EQR > 1, 1, 
                               EQR)))                 %>%
    select(-TYPO_I2M2, -RESULTAT, -WORST, -BEST)      %>%
    mutate(METRIQUE = factor(METRIQUE, 
                             levels = c("SHANNON", "ASPT", "POLYVOLTIN",
                                        "OVOVIVIPARE", "RICHESSE")
    ))                       %>%
    spread(key = METRIQUE, value = EQR)
}

## Fonction permettant le calcul de l'I2M2
funI2M2 <- function(Table, De) {
  gather(Table, key = METRIQUE, value = EQR, 
         -CODE_OPERATION)                     %>%
    left_join(gather(De, key = PRESSION, value = DE, 
                     -METRIQUE),
              by = "METRIQUE")                %>%
    group_by(CODE_OPERATION, PRESSION)        %>%
    summarise(i2m2 = sum(EQR * DE) / sum(DE)) %>%
    group_by(CODE_OPERATION)                  %>%
    summarise(I2M2 = mean(i2m2))
}

## Fonction permettant de formatter les resultats pour le fichier de sortie
funFormat <- function(Table) {
  codes <- structure(c(8058, 8057, 8056, 8055, 8054, 7613, 8050), 
                     names = c("SHANNON", "ASPT", "POLYVOLTIN", 
                               "OVOVIVIPARE", "RICHESSE", "I2M2",
                               "nbContributif"))
  
  gather(Table, key = METRIQUE, value = RESULTAT, 
         -CODE_OPERATION)                   %>%
    mutate(CODE_PAR = codes[METRIQUE])      %>%
    filter(!is.na(CODE_PAR))                %>%
    mutate(RESULTAT = as.numeric(RESULTAT)) %>%
    select(CODE_OPERATION, CODE_PAR, RESULTAT)
}


## Fonction generant les commentaires
funCommentaires <- function(Table, TaxaContributifs) {
  bind_rows(
    group_by(Table, CODE_OPERATION) %>%
      summarise(COMMENTAIRE = 
                  ifelse(RESULTAT[CODE_PAR == 8050] > 0 &
                           is.na(RESULTAT[CODE_PAR == 8058]),
                         "Type de cours d'eau non pris en compte",
                         "")) %>%
      mutate(CODE_PAR = 8058) %>%
      filter(COMMENTAIRE != ""),
    group_by(TaxaContributifs, CODE_OPERATION) %>%
      summarise(COMMENTAIRE = 
                  ifelse(nbNonContributif > 0,
                         paste0("Les taxons suivants, representant ",
                                round(100 * nbNonContributif / nbTotal),
                                "% des taxons du prelevement, n'ont pas",
                                " ete pris en compte dans le calcul : ",
                                taxaNonContributifs, "."),
                         "")) %>%
      mutate(CODE_PAR = 8050) %>%
      filter(COMMENTAIRE != "")
  )
}

## Fonction permettant de faire les arrondis a l'inferieur si 0 a 4 et au superieur si 5 a 9
funArrondi <- function (x, digits = 0) {
  .local <- function(x, digits) {
    x <- x * (10^digits)
    ifelse(abs(x%%1 - 0.5) < .Machine$double.eps^0.5,
           ceiling(x)/(10^digits),
           round(x)/(10^digits))
  }

  if (is.data.frame(x))
    return(data.frame(lapply(x, .local, digits)))
  .local(x, digits)
}

## Fonction initialisant le fichier de sortie
funSortie <- function(data_entree, paramsOut, ...) {
  select(data_entree, ...) %>%
    distinct()             %>%
    (function(df) {
      df[rep(1:nrow(df), each = nrow(paramsOut)),] %>%
        as.tbl()
    })                     %>%
    mutate(CODE_PAR = rep(paramsOut$CODE_PAR,
                          n() / nrow(paramsOut)),
           LIB_PAR  = rep(paramsOut$LIB_PAR,
                          n() / nrow(paramsOut)))
}

## Fonction permettant d'ecrire le fichier de sortie
funResult 		<- function(indic, vIndic, heure_debut,
                        data_sortie, data_complementaire, complementaire,
                        file, file_complementaire)
{
    # determination du temps de calcul
    heure_fin       <- Sys.time()
    heure_dif       <- heure_fin - heure_debut
    temps_execution <- paste0(round(heure_dif, 2),
                              attr(heure_dif, "units"))

    # creation du bandeau d'information
    etiquette <- paste(indic, vIndic, Sys.Date(),
                       "Temps d'execution :", temps_execution,
                       sep = ";")

    # sortie du bandeau d'information
    cat(paste0(etiquette, "\n"), file = file, sep = "")

    # sortie du fichier de sortie
    write.table(data_sortie, row.names = FALSE, quote = FALSE, sep = ";",
                file = file, append = TRUE)

    # Sortie complementaire
    if(complementaire)
    {
        if (file == "") {
            print("Fichier")
        }

        cat(paste0(etiquette, "\n"), file = file_complementaire, sep = "")
        write.table(data_complementaire, row.names = FALSE, quote = FALSE,
                    sep = ";", file = file_complementaire, append = TRUE)
    }

}# fin de la fonction funResult

## INITIALISATION DU TRAITEMENT ----
# Ne pas afficher les messages d'avis ni d'erreur
options(warn = -1)

File           <- "I2M2_entree_op100.txt"
complementaire <- TRUE

# Initialisation de l'heure
heure_debut <- Sys.time()

##  IMPORT DES FICHIERS ----
# Import du fichier d'entree
data_entree <- read.table(File, header = TRUE, sep = "\t",
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

fichierResultat               <- paste0(indic, "_", vIndic, "_resultats.csv")
fichierResultatComplementaire <- paste0(indic, "_", vIndic,
                                        "_resultats_complementaires.csv")
funResult(indic               = indic,
          vIndic              = vIndic,
          heure_debut         = heure_debut,
          data_sortie         = data_sortie,
          data_complementaire = data_complementaire,
          complementaire      = complementaire,
          file                = fichierResultat,
          file_complementaire = fichierResultatComplementaire)
