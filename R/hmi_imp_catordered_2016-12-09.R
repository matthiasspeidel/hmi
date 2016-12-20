####################################################################
###Funktion zur Imputation/Synthetisierung kategorialer Variablen###
####################################################################

#gehe davon aus, dass Y kategorial ist.
#unter der Annahme, dass sich die Wkeit ob eine Beobachtung
#in Kategorie a, b, c  oder sonst eine kommt, von X abh?ngig ist,
#k?nnte ich eigentlich eine http://de.wikipedia.org/wiki/
#Multinomiale_logistische_Regression
#sch?tzen... Das Problem ist, dass f?r Variablen mit vielen Kategorien
#die einzelnen Kategorien evtl sehr schwach besetzt sind und damit das
#Separationsproblem auftritt. Deshalb CART-Ansatz


imp_orderedcat <- function(data.imp.ordered.categorical, y.variable.name, M){

  y.org <- subset(data.imp.ordered.categorical, select = y.variable.name)

  Kategorien <- levels(data.imp.ordered.categorical[, y.variable.name]) #as.numeric()
  #nach M?glichkeit, sollte daten ein data.frame sein

  #Initialisierung des Vektors mit den (zuk?nftig) imputierten
  #bzw. synthetisierten Werten.
  y.imp <- y.org

  #the missing indactor indicates, which values of y are missing.
  mis.indicator <- is.na(y.org)
  X.org <- subset(data.imp.ordered.categorical, select = -which(names(data.imp.ordered.categorical) == y.variable.name))

  #MS: schaue welche Variablen ein strukturelles/Filterbedingtes NA haben. Bsp: Alter vom Kind = NA, wenn Anzahlkinder = 0.
  mis.indicator.X.org <- apply(X.org, 2, FUN = function(x) any(is.na(x)))

  #MS: entferne diese Variablen
  X.org <- subset(X.org, select = !mis.indicator.X.org)

  #ermittle die Variablen, die Faktorvariablen sind
  #!!!stand heute gehe ich von nur einer Faktorvariablen aus!!!!
  factor.variables = unlist(lapply(data.imp.ordered.categorical, is.factor))


  level.lengths = unlist(lapply(data.imp.ordered.categorical, function(x) length(levels(x))))

  problematic.factor.variables = level.lengths > 25 & names(data.imp.ordered.categorical) != y.variable.name

  #alternative declaration in one step:
  #problematic.factor.variables =
  #  unlist(lapply(data.imp.ordered.categorical[, factor.variables],
  #                function(x) length(levels(x)))) > 25 & names(data.imp.ordered.categorical[, factor.variables]) != y.variable.name


  #get an R-variable, that indicates the length of that problematic factor-variable
  problematic.factor.variables.length <- level.lengths[problematic.factor.variables]

  #check if there is an problem
  is.problem = any(problematic.factor.variables)

  #How many steps do I need, to get clusters with less than 25 factorlevels?

  #in the case, that no problem exists (a factor variable has a problematic length)...
  if(!is.problem){
    #... for the later use it is set to 1
    problematic.factor.variables.length = 1
  }

  #in wieviele Gruppen muss ich die Faktorvariable aufteilen,
  #um gruppen <= 25 zu kriegen?
  #ceiling(problematic.factor.variables.length/25)
  nsteps <- ceiling(problematic.factor.variables.length/25)
  #wenn es insgesamt weniger als 25 factorstufen gibt, sollte das 1 sein,
  #und ganz normal funktionieren.


  for(l1 in 1:nsteps){


    if(is.problem){
      #f?r gleich gro?e Gruppen f?nde ich es gut, wenn die Aufteilung
      current.levels = (floor(problematic.factor.variables.length * (l1-1)/nsteps) + 1):
        floor(problematic.factor.variables.length * l1/nsteps)


      #nimm nur die Beobachtungen die in den ersten 25 Faktorstufen vorkokmmen
      current = data.imp.ordered.categorical[, problematic.factor.variables] %in%
        levels(data.imp.ordered.categorical[, problematic.factor.variables])[current.levels]
      #in dem unproblematischen Fall, ist current dann einfach ein Vektor mit
      #nrow(data.imp.ordered.categorical) TRUEs. Also:
      #current <- array(TRUE, dim = nrow(data.imp.ordered.categorical))

      #data.imp.ordered.categorical[current,] kann noch nicht verwendet werden, da die Faktorvariable
      #noch nicht upgedated ist (und R nicht wei?, dass nun nur noch 25 Stufen
      #vertreten sind)
      data.current = data.imp.ordered.categorical[current, ]
      data.current[, problematic.factor.variables] =
        as.factor(as.character(data.current[, problematic.factor.variables]))
    }else{
      current <- TRUE
      data.current <- data.imp.ordered.categorical[current, ]
    }

    ry <- !is.na(y.org)
    #MS: VERSUCH: steige fr?her bei mice in
    tmpmice.data <- cbind(y.org, X.org)
  everything <- mice(data = tmpmice.data, m = M,
              method = "polr",
              predictorMatrix = (1 - diag(1, ncol(tmpmice.data))),
              visitSequence = (1:ncol(tmpmice.data))[apply(is.na(tmpmice.data),2,any)],
              post = vector("character", length = ncol(tmpmice.data)),
              defaultMethod = "polr",
              maxit = 10,
              diagnostics = TRUE,
              printFlag = FALSE,
              seed = NA,
              imputationMethod = NULL,
              defaultImputationMethod = NULL,
              data.init = NULL)

    y.imp <- y.org
    y.imp[!ry] <- everything$imp[[1]][, 1]

    #MS: in mice.impute.polr() macht die Funktion augment() Probleme mit dem Befehl  sort(unique(unclass(y)))
    #MS: Der (so verstehe ich (MS) es heute (14.01.2015)) bestimmen soll, wieviele Kategorien in der Variable enthalten sind.
    #MS: denn die L?nge des Objekts mit mit einer Maximalanzahl verglichen und bei ?berschreitung die Imputation abgebrochen.


    }
  ret <- y.imp
  return(ret)
}




