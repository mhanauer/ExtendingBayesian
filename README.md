# Jags-Ymet-XmetSsubj-MrobustHierQuadWtMatt'sTest.R 
# Accompanies the book:
#  Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: 
#  A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier.

source("DBDA2E-utilities.R")

#===============================================================================

genMCMC = function( data , xName="x" , yName="y" , sName="s" , wName=NULL ,
                    numSavedSteps=10000 , thinSteps = 1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault ) { 

  #-----------------------------------------------------------------------------
  # THE DATA.
  y = data[,yName]
  x = data[,xName]
  s = as.numeric(data[,sName])
  if ( !is.null(wName) ) {
    w = data[,wName]
  } else {
    w = rep(1,length(y))
  }
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  #Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    s = s ,
    w = w ,
    Nsubj = max(s)  # should equal length(unique(s))
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  # Standardize the data:
  data {
    Ntotal <- length(y)
    xm <- mean(x)
    ym <- mean(y)
    wm <- mean(w)
    xsd <- sd(x)
    ysd <- sd(y)
    for ( i in 1:length(y) ) {
      zx[i] <- ( x[i] - xm ) / xsd
      zy[i] <- ( y[i] - ym ) / ysd
      zw[i] <- w[i] / wm 
    }
  }
  # Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      zy[i] ~ dt( zbeta0[s[i]] + zbeta1[s[i]] * zx[i] + zbeta2[s[i]] * zx[i]^2 , 
                  1/(zw[i]*zsigma)^2 , nu )
    }
    for ( j in 1:Nsubj ) {
      zbeta0[j] ~ dnorm( zbeta0mu , 1/(zbeta0sigma)^2 )  
      zbeta1[j] ~ dnorm( zbeta1mu , 1/(zbeta1sigma)^2 )
      zbeta2[j] ~ dnorm( zbeta2mu , 1/(zbeta2sigma)^2 )
    }
    # Priors vague on standardized scale:
    zbeta0mu ~ dnorm( 0 , 1/(10)^2 )
    zbeta1mu ~ dnorm( 0 , 1/(10)^2 )
    zbeta2mu ~ dnorm( 0 , 1/(10)^2 )
    zsigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbeta0sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbeta1sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    zbeta2sigma ~ dunif( 1.0E-3 , 1.0E+3 )
    nu ~ dexp(1/30.0)
    # Transform to original scale:
    for ( j in 1:Nsubj ) {
      beta2[j] <- zbeta2[j]*ysd/xsd^2
      beta1[j] <- zbeta1[j]*ysd/xsd - 2*zbeta2[j]*xm*ysd/xsd^2  
      beta0[j] <- zbeta0[j]*ysd  + ym - zbeta1[j]*xm*ysd/xsd + zbeta2[j]*xm^2*ysd/xsd^2 
    }
    beta2mu <- zbeta2mu*ysd/xsd^2
    beta1mu <- zbeta1mu*ysd/xsd - 2*zbeta2mu*xm*ysd/xsd^2  
    beta0mu <- zbeta0mu*ysd  + ym - zbeta1mu*xm*ysd/xsd + zbeta2mu*xm^2*ysd/xsd^2 
    sigma <- zsigma * ysd
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta1" ,  "beta2" ,
                  "beta0mu" , "beta1mu" , "beta2mu" ,
                  "zbeta0" , "zbeta1" , "zbeta2" ,
                  "zbeta0mu" , "zbeta1mu" , "zbeta2mu" ,
                  "sigma" , "nu" , 
                  "zsigma", "zbeta0sigma" , "zbeta1sigma", "zbeta2sigma" )
  adaptSteps = 1000  # Number of steps to "tune" the samplers
  burnInSteps = 10000 
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=adaptSteps ,
                          burnin=burnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]

  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#===============================================================================

smryMCMC = function(  codaSamples , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=FALSE)
  paramNames = colnames(mcmcMat)
  summaryInfo = NULL
  for ( pName in paramNames ) {
    summaryInfo = rbind( summaryInfo ,  summarizePost( mcmcMat[,pName] ) )
  }
  rownames(summaryInfo) = paramNames
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}
