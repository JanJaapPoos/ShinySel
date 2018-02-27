library(shiny)
library(deSolve)

ui <- fluidPage(titlePanel("Selectivity example in shiny"),
sidebarLayout(
  sidebarPanel("inputs",
               numericInput("meshsize", "meshsize", min = 1, max = 150,step=1,   value = c(80)),
               textInput("selple", "Selectivity plaice (sf=2.24,sr=3.66)", 
                            value="2.24,3.66"),
               textInput("vblgple", "VBLG Growth plaice (Linf=50,k=0.3)",
                           value = "50,0.3"),
               textInput("selsol", "Selectivity sole (sf=3.34,sr=3.55)", 
                           value = "3.34,3.55"),
               textInput("vblgsol", "VBLG Growth sole (Linf=30,k=0.2)",
                           value = "30,0.2")),
  mainPanel( plotOutput("coolplot", height="750px"),verbatimTextOutput("results"))
))

server <- function(input, output) {
  
  
  output$results <- renderText({ 
    "test but should be numerical outcomes of YPR curve"
    })
  
  
  output$coolplot <- renderPlot({
  
    nselple <- as.numeric(unlist(strsplit(input$selple,",")))
    nselsol <- as.numeric(unlist(strsplit(input$selsol,",")))
    nvblgple <- as.numeric(unlist(strsplit(input$vblgple,",")))
    nvblgsol <- as.numeric(unlist(strsplit(input$vblgsol,",")))
    
    sels <- function (meshsize, SF,SR, prange,MLS=0) {
      L50 <- SF*(meshsize/10) #division by 10 is needed for mm to cm
      res <- plogis(prange, L50,SR) 
      res[(prange < MLS)] <- 0   
      return(res)
    }
    
    age <- seq(0,15,0.01)
    alphaple <- 0.001
    alphasol <- 0.001
    betaple <- 3
    betasol <- 3
    prange <- seq(0, 50, by = 0.01)
  
    
     outple <- sels(input$meshsize,nselple[1],nselple[2],prange)
     outsol <- sels(input$meshsize,nselsol[1],nselsol[2],prange)
     
     par(mfrow=c(5,2), mar=c(4,4,0.1,0.1), mgp=c(2,1,0))
     plot(x=prange, y=outple, type = "l", xlab = "Fish length (cm)", ylab = "population", ylim=c(0,max(c(outple,outsol))), first.panel=grid())
     text(x=5,y=0.9,"ple")
     abline(v=nselple[1]*(input$meshsize/10))
     abline(v=27, lty=2, col="red")
     lines(x=prange[outple>0.25 &outple <0.75], y=rep(0.01,length(prange[outple>0.25 &outple <0.75])))
     text(x=prange[outple>0.25 &outple <0.75][1],y=0.05,rev(prange[outple>0.25 &outple <0.75])[1]-prange[outple>0.25 &outple <0.75][1])
   
     plot(x=prange, y=outsol, type = "l", xlab = "Fish length (cm)", ylab = "population", ylim=c(0,max(c(outple,outsol))), first.panel=grid())
     text(x=5,y=0.9,"sol")
     abline(v=nselsol[1]*(input$meshsize/10))
     abline(v=24, lty=2, col="red")
     lines(x=prange[outsol>0.25 &outsol <0.75], y=rep(0.01, length(prange[outsol>0.25 &outsol <0.75])))
     text(x=prange[outsol>0.25 &outsol <0.75][1],y=0.05,rev(prange[outsol>0.25 &outsol <0.75])[1]-prange[outsol>0.25 &outsol <0.75][1])
     
     
     lenple <- nvblgple[1] * (1-exp(- nvblgple[2] * age))
     lensol <- nvblgsol[1] * (1-exp(- nvblgsol[2] * age))
     
     wtsple <- alphaple*(lenple)^betaple
     wtssol <- alphasol*(lensol)^betasol
     
     plot(x=age,y=wtsple, type="l", ylim=c(0,max(c(wtsple,wtssol))))
     text(x=5,y=0.9,"ple")
      
     plot(x=age,y=wtssol, type="l")
     text(x=5,y=0.9,"sol")
     
     ageselple <- sels(input$meshsize,nselple[1],nselple[2],lenple)
     ageselsol <- sels(input$meshsize,nselsol[1],nselsol[2],lensol)
     
     ageselpleland <- sels(input$meshsize,nselple[1],nselple[2],lenple, MLS=27)
     ageselsolland <- sels(input$meshsize,nselsol[1],nselsol[2],lensol, MLS=24)
     
     plot(x=age, y=ageselple, type="l", ylim=c(0,1))
     lines(x=age, y=ageselpleland, lty=2)
     
     plot(x=age, y=ageselsol, type="l", ylim=c(0,1))
     lines(x=age, y=ageselsolland, lty=2)
     
     mat <- rep(1, length(age))
     M   <- rep(0.001, length(age))
     
     # The procedure for calculating yield curve is similar to yield curve. 
     # However, we have to account for effect of S-R curve
     # We make use of: SSB, = Spawner biomass per recruit (SSBPR)  * recruitment (R): (SSB=SSPR * R)
     # next we can substitute this in Beverton and Holt curve: SSB/(alpha + beta *SSB)
     #
     # R=  (SSBPR *R)/ (alpha + beta *SSBPR *R) = (SSPR- alpha)/(beta*SSPR)
     #
     # Then 1) Calculate the YPR and SSBPR as a function of exploitation rate
     #      2) Calculate R given SSBPR and the stock-recruitment relationship.
     #      3) Multiply YPR and SSBPR by R to calculate yield and spawner biomass.
     
     # first we need SR parameters for BH in R= (SSBPR *R)/ (alpha + beta *SSBPR *R) formulation
     
     # SSB_SR <- SSB[-length(SSB)]
     #  R_SR   <- recruitment[-1]
     
     # logAlpha  <- log(0.5)
     #  logBeta   <- log(0.001)
     # logSigma <- 6.5 
     
     # par <- c(logAlpha,logBeta,logSigma)
     
     #  logLikelihood <- function(par){
     #   Rexpected <- SSB_SR/(exp(par[1]) + exp(par[2]) * SSB_SR ) 
     #    nll <- -sum(dnorm(R_SR,Rexpected,exp(par[3]),TRUE))
     #    return(nll)
     #  } 
    
     #res contains all important information from optimizer, including MLE estimates
     # res <- optim(par, logLikelihood, method = "BFGS", control=list(trace=T, maxit=500))
      
     # alpha <- exp(res$par[1])
     #  beta  <- exp(res$par[2])   
     
     #  plot(x=SSB_SR,y=R_SR, xlim=c(0,2500), ylim=c(0,2500))
     #  lines(x=0:2500, y=(0:2500)/(alpha + beta * 0:2500 ))
      
     # Next we will perform YPR curve for ple described at the top of the page
     vals  <- seq(0,0.02,0.0005)
     CPRple <- LPRple <- SSBPRple <- eqYple <- eqSSBple <- Fbarple <- numeric(length(vals))
     ix <- 1
     for (ii in vals){
       Z         <- M + (ageselple * ii)
       N         <- c(1, exp(-cumsum(Z)))[1:length(age)]
       Catch     <- ((ageselple * ii )/Z) * N * (1-exp(-Z))
       Land      <- ((ageselpleland * ii )/Z) * N * (1-exp(-Z))
       CPRple[ix]   <- sum(Catch * wtsple)             # Yield per recruit
       LPRple[ix]   <- sum(Land * wtsple)             # Yield per recruit
       SSBPRple[ix] <- sum(N * wtsple *mat)     # SSBPR
    #  R         <- (SSBPR[ix]- alpha)/(beta*SSBPR[ix]) # equilibrium R
    #  eqSSB[ix] <- SSBPR[ix] * R                       # equilibrium SSB in curve given F  
    #   eqY[ix]   <- YPR[ix] * R                         # equilibrium Yield in curve given F
        Fbarple[ix]  <- mean((ageselple * ii))
        ix <- ix +1
     } 
     
     CPRsol <- LPRsol <- SSBPRsol <- eqYsol <- eqSSBsol <- Fbarsol <- numeric(length(vals))
     ix <- 1
     for (ii in vals){
       Z         <- M + (ageselsol * ii)
       N         <- c(1, exp(-cumsum(Z)))[1:length(age)]
       Catch     <- ((ageselsol * ii )/Z) * N * (1-exp(-Z))
       Land      <- ((ageselsolland * ii )/Z) * N * (1-exp(-Z))
       CPRsol[ix]   <- sum(Catch * wtssol)             # Yield per recruit
       LPRsol[ix]   <- sum(Land * wtssol)             # Yield per recruit
       SSBPRsol[ix] <- sum(N * wtssol *mat)     # SSBPR
       #    R         <- (SSBPR[ix]- alpha)/(beta*SSBPR[ix]) # equilibrium R
       #    eqSSB[ix] <- SSBPR[ix] * R                       # equilibrium SSB in curve given F  
       #    eqY[ix]   <- YPR[ix] * R                         # equilibrium Yield in curve given F
       Fbarsol[ix]  <- mean((ageselsol * ii))
       ix <- ix +1
     } 
     
  #   ########################### where is FMSY?
  # #  FMSY_BH <- Fbar[which(eqY == max(eqY))]
  #   
  #   # make plots where we can compare YPR en Yield curve, and SSBPR and eqSSB curve
     plot(Fbarple,CPRple, type="l")
     lines(Fbarple,LPRple, type="l")
     plot(Fbarsol,CPRsol, type="l")
     lines(Fbarsol,LPRsol, type="l")
  # #  plot(Fbar,eqY, type="l")
  # #  abline(v=FMSY_BH, lty=2, col="magenta")
  # #  text(1,300,paste("BH FMSY",round(FMSY_BH,3)))
  #   
  # #  plot(Fbar,SSBPR, type="l")
  # #  plot(Fbar,eqSSB, type="l")
  # #  abline(v=FMSY_BH, lty=2, col="magenta")
  #   
  }, res=100)
}
shinyApp(ui = ui, server = server)
