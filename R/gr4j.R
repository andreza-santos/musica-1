#' Estimation of potential evapotranspiration based on temperature
#'
#' @param TAS Temperature
#' @param DTM datum
#' @param lat latitude
#'
#' @return
#' @export
#'
#' @examples
t.pet=function(TAS, DTM, lat = 49.8){

  Ra=function(lat, J){

    rlat = pi/180 * lat
    Gsc = 0.0820
    dr = 1 + 0.033 *	cos(J*2*pi/365)
    delta = 0.409 * sin(J*2*pi/365-1.39)
    om = acos (-tan(rlat) * tan(delta))
    (24*60)/pi * Gsc * dr * (om * sin(rlat)*sin(delta) + cos(rlat)*cos(delta)*sin(om))
  }

  J = as.integer(format(as.Date(DTM),'%j'))
  Re = Ra(rep(lat, length(J)), J)
  t5 = TAS+5
  #t5[t5<=5]=0

  pet = 0.408 * Re * t5 / 100
  pet[pet<0] = 0
  pet
}


#' Runoff simulation with the GR4J model
#'
#' @param data data.table with columns DTM (datum), PR (precipitation), TAS (temperature)
#' @param x parameters of the GR4J model
#' @param lat latitude
#' @return estimated runoff
#' @export
#'
#' @examples
#' @useDynLib musica sma_gr4j
#' @useDynLib musica routing_gr4j
gr4j = function(data, x = c(350, 0, 40, 0.5), lat = 49.8, return_components = FALSE){

  dta = copy(data)
  dta[, E:=t.pet(TAS, DTM, lat = lat)]
  er=gr4j.sim(cbind(P=dta$PR, E=dta$E),x1=x[1],etmult = 1, S_0 = 0.5, return_state = FALSE)
  r = gr4jrouting.sim(er,x[2],x[3],x[4], eps=0.001, R_0 = .5, return_components = FALSE)
  r[1] = r[2] # to replace NA at the first time step (inherent to GR4J)
  r
}

gr4j.sim <-  function(DATA,
                      x1, etmult = 1, S_0 = 0.5,
                      return_state = FALSE)
{
  stopifnot(c("P","E") %in% colnames(DATA))
  ## check values
  stopifnot(x1 >= 0)
  stopifnot(etmult >= 0)
  stopifnot(S_0 >= 0)
  stopifnot(S_0 <= 1)

  inAttr <- attributes(DATA[,1])
  DATA <- as.ts(DATA)
  P <- DATA[,"P"]
  E <- DATA[,"E"] * etmult

  ## skip over missing values (maintaining the state S)
  bad <- is.na(P) | is.na(E)
  P[bad] <- 0
  E[bad] <- 0
  #   COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
  #   if (COMPILED) {
  ans <- .C('sma_gr4j',
            as.double(P),
            as.double(E),
            as.integer(NROW(DATA)),
            as.double(x1),
            as.double(S_0),
            U = double(NROW(DATA)),
            S = double(NROW(DATA)),
            ET = double(NROW(DATA)),
            NAOK=FALSE, DUP=FALSE)
  U <- ans$U
  S <- ans$S
  ET <- ans$ET
  attributes(U) <- inAttr
  ## re-insert missing values
  U[bad] <- NA
  ans <- U
  if (return_state) {
    attributes(S) <- attributes(ET) <- attributes(U)
    ans <- cbind(U=U, S=S, ET=ET)
  }
  return(ans)
}

gr4jrouting.sim <-
  function(U, x2, x3, x4, R_0 = 0,
           return_components = FALSE,
           epsilon = hydromad.getOption("sim.epsilon"))
  {
    ## check values
    stopifnot(is.numeric(x2))
    stopifnot(x3 >= 0)
    stopifnot(x4 >= 0.5)
    stopifnot(R_0 >= 0)
    stopifnot(R_0 <= 1)

    inAttr <- attributes(U)
    U <- as.ts(U)

    n <- ceiling(x4)
    m <- ceiling(x4 * 2)

    ## S-curves: cumulative proportion of input with time
    SH1 <- (1:x4 / x4) ^ (5/2)
    SH2 <- c(0.5 * SH1,
             1 - 0.5 * (2 - floor(x4+1):(2*x4) / x4) ^ (5/2)
    )
    ## unit hydrographs
    UH1 <- diff(SH1)
    if (length(UH1)==0) UH1=1
    UH2 <- diff(SH2)

    ## skip over missing values (maintaining the state)
    bad <- is.na(U)
    U[bad] <- 0

    Q9 <- filter(0.9 * U, UH1, sides = 1)
    Q1 <- filter(0.1 * U, UH2, sides = 1)

    ## fill in values undefined by UH filter
    Q9[seq_along(UH1)] <- 0
    Q1[seq_along(UH2)] <- 0
    bad[seq_along(UH2)] <- TRUE

    #   COMPILED <- (hydromad.getOption("pure.R.code") == FALSE)
    #   if (COMPILED) {
    ans <- .C('routing_gr4j',
              as.double(Q9),
              as.double(Q1),
              as.integer(length(U)),
              as.double(x2),
              as.double(x3),
              as.double(R_0),
              Qr = double(length(U)),
              Qd = double(length(U)),
              R = double(length(U)),
              NAOK=FALSE, DUP=FALSE)
    Qr <- ans$Qr
    Qd <- ans$Qd
    R <- ans$R

    ## zap simulated values smaller than epsilon
    Qr[Qr < epsilon] <- 0
    Qd[Qd < epsilon] <- 0

    attributes(Qr) <- attributes(Qd) <- attributes(R) <- inAttr

    ## re-insert missing values
    Qr[bad] <- NA
    Qd[bad] <- NA

    if (return_components) {
      return(cbind(Xr = Qr, Xd = Qd, R = R))
    } else {
      return(Qr + Qd)
    }
  }

