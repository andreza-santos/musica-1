
ext_cast = function(out, what){
  dcast.data.table(out$dta[RID==what, .(variable, DTM, D1)], DTM ~ variable, value.var = 'D1')
}


approx = function(x, y = NULL, xout, method = "linear", n = 50, yleft, yright, rule = 1, f = 0, ties = mean){

  if (length(unique(x)) == 1)  {
    return(list(x = xout, y = rep(mean(y), length(xout))))
  } else {
    return(
      stats::approx(x, y, xout, method , n , yleft, yright, rule, f, ties)
    )
  }
}

msTrans_anom = function(dta, agg_by = month, wet_int_thr = 0.1, maxiter = 10, tol = 1e-4, report = FALSE, ipreci = .001, year_starts = months(0), scales = c('G1', 'Y1', 'M3', 'M1', 'D1')){

  giveAnom = Vectorize(function(i){

    mdta[, cut := as.IDate(musica:::cut.Date(DTM, breaks = cpr[i])), by = .(RID, variable)]
    mdta[, GRP := month(cut %m-% year_starts)]

    if (i == 1){
      mdta[, G1 := mean(ORIG), by = .(RID, variable)]
    } else {
      tex = paste0('c(', paste0(pr[1:(i-1)], '[1]', collapse = ',') , ')')

      mdta[, pr[i] := dif(mean(ORIG), rev_difv(eval(parse(text = tex )), variable[1]), variable[1]), by = .(RID, variable, cut)]
    }
    return(NULL)
  })


  pr = scales
  cpr = code2period(pr)
  ip = seq(0, 1, by = ipreci)

  pdif = c(NA, NA)
  cdif = c(NA, NA)
  pp = c(NA, NA)
  idif = pdif - cdif
  ite = 1

  mdta = melt(dta, id.vars = c('RID', 'DTM'), value.name = 'ORIG')


  while ( (ite <= maxiter) & ((max(abs(idif)) > tol) | any(is.na(idif))) ) {

  cat('\n', format(ite, width =3), '(', format(signif(idif, 5), width = 12), format(pp, width = 6), ') :', '\t')

    invisible(giveAnom(1:length(pr)))

    for (i in 1:length(pr)){

      invisible(giveAnom(i))

      mdta[, MEA := eval(parse(text = pr[i]))]

      smry = mdta[, .(MEA = MEA[1]), by = .(RID, variable, cut, GRP)]

      f = list()
      for (g in smry[, unique(GRP)]){
        for (v in smry[, unique(variable)]){
          f[[length(f)+1]] = #if (cpr[i]!=0) {
            smry[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(fitQmap(obs = MEA[RID=='TO'], mod = MEA[RID=='FROM'], method = 'QUANT', wet.day = FALSE, qstep = .001))) ]
            #cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(approxfun(x = p, y = CF, rule = 2)))]
         # } else {
          #  cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(apxfun( CF )))]
          #}

        }
      }

      f = rbindlist(f)

      mdta = f[mdta, on = c('GRP', 'variable')]
      mdta[RID=='MODIF', pr[i] := doQmapQUANT.default( eval(parse(text = pr[i])), fun[[1]] ), by = .(GRP, variable) ]
      mdta[RID=='FROM', pr[i] := doQmapQUANT.default( eval(parse(text = pr[i])), fun[[1]] ), by = .(GRP, variable) ]

      mdta[, c('GRP', 'fun', 'cut', 'MEA') := NULL]
      tex = paste0('c(', paste0(scales, '[1]', collapse = ',') , ')')
      mdta[, ORIG := rev_difv(eval(parse(text = tex)), variable[1]), by = .(RID, variable, DTM)]

            #
      # #
      # mdta[, mean(ORIG), by = .(variable, RID)]
      # mdta[, mean(ORIG), by = .(variable, RID, year(DTM))][, sd(V1), by = .(variable, RID)]
      #
      # mo = mdta[, mean(ORIG), by = .(variable, RID, year(DTM), month(DTM))]
      # ggplot(mo[RID %in% c('FROM', 'TO'), .(V1 = sort(V1), year = 1:length(V1)), by = .(month, variable, RID)]) +  geom_line(aes(x = year, y = V1, col = RID))+  geom_point(aes(x = year, y = V1, col = RID)) + facet_grid(variable ~ month, scale = 'free')
      # #a = mdta[, .(ks.test(G1[RID == 'FROM'], G1[RID == 'TO'])$p.value, ) ]


      }

    invisible(giveAnom(1:length(pr)))

    err = data.table()
    for (i in 1:length(pr)){
      suppressWarnings({
      err = rbind(err, mdta[, .(id = c('D', 'p.value', 'alternative', 'method', 'call'), scale = pr[i], unlist(ks.test(eval(parse(text = pr[i]))[RID == 'FROM'], eval(parse(text = pr[i]))[RID == 'TO']) )), by = variable ])
    })}

    cerr = dcast(err, variable + scale ~id, value.var = 'V3')
    DD = cerr[, sum(as.double(D)), by = variable][, V1]
    pp = cerr[, sum(as.numeric(p.value) > .1), by = variable][, V1]

    idif = pdif - DD
    pdif = DD
    ite = ite + 1

    }

  cat('\n')
  list(dta = copy(mdta[, .(variable, RID, DTM, D1 = ORIG)]))


}


######
msTrans_abs = function(dta, agg_by = month, wet_int_thr = 0.1, maxiter = 10, tol = 1e-4, report = FALSE, ipreci = .001, year_starts = months(0), scales = c('G1', 'Y1', 'M3', 'M1', 'D1')){

  giveAbs = Vectorize(function(i){

    mdta[, cut := as.IDate(musica:::cut.Date(DTM, breaks = cpr[i])), by = .(RID, variable)]
    mdta[, GRP := month(cut %m-% year_starts)]
    mdta[, pr[i] := mean(ORIG), by = .(RID, variable, GRP, cut)]
    return(NULL)
  })



  apxfun = function(y){
    function(v) y
  }

  pr = scales
  cpr = code2period(pr)

  pdif = c(NA, NA)
  cdif = c(NA, NA)
  pp = c(NA, NA)
  idif = pdif - cdif
  ite = 1

  mdta = melt(dta, id.vars = c('RID', 'DTM'), value.name = 'ORIG')


  while ( (ite <= maxiter) & ((max(abs(idif)) > tol) | any(is.na(idif))) ) {

    cat('\n', format(ite, width =3), '(', format(signif(idif, 5), width = 12), format(pp, width = 6), ') :', '\t')
    invisible(giveAbs(1:length(pr)))

    for (i in 1:length(pr)){

      invisible(giveAbs(i))

      mdta[, MEA := eval(parse(text = pr[i]))]

      smry = mdta[, .(MEA = MEA[1]), by = .(RID, variable, cut, GRP)]

      f = list()
      for (g in smry[, unique(GRP)]){
        for (v in smry[, unique(variable)]){
          f[[length(f)+1]] = #if (cpr[i]!=0) {
            smry[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(fitQmap(obs = MEA[RID=='TO'], mod = MEA[RID=='FROM'], method = 'QUANT', wet.day = FALSE, qstep = .001))) ]
          #cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(approxfun(x = p, y = CF, rule = 2)))]
          # } else {
          #  cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(apxfun( CF )))]
          #}

        }
      }

      f = rbindlist(f)

      mdta = f[mdta, on = c('GRP', 'variable')]

      mdta[RID %in% c('MODIF', 'FROM'), CORR := doQmapQUANT.default( eval(parse(text = pr[i])), fun[[1]] ), by = .(GRP, variable, RID) ]

      mdta[RID %in% c('MODIF', 'FROM'), CF := dif(CORR, eval(parse(text = pr[i])), variable[1]), by = .(GRP, variable, RID) ]
      mdta[RID %in% c('MODIF', 'FROM'), ORIG:= rev_dif(CF, ORIG, variable[1]), by = .(GRP, variable, RID)]
      mdta[, c('GRP', 'fun', 'CORR', 'CF', 'MEA', 'cut'):=NULL]
      #
      # #
      # mdta[, mean(ORIG), by = .(variable, RID)]
      # mdta[, mean(ORIG), by = .(variable, RID, year(DTM))][, sd(V1), by = .(variable, RID)]
      #
      # mo = mdta[, mean(ORIG), by = .(variable, RID, year(DTM), month(DTM))]
      # ggplot(mo[RID %in% c('FROM', 'TO'), .(V1 = sort(V1), year = 1:length(V1)), by = .(month, variable, RID)]) +  geom_line(aes(x = year, y = V1, col = RID))+  geom_point(aes(x = year, y = V1, col = RID)) + facet_grid(variable ~ month, scale = 'free')
      # #a = mdta[, .(ks.test(G1[RID == 'FROM'], G1[RID == 'TO'])$p.value, ) ]


    }

    giveAbs(1:length(pr))

    err = data.table()
    for (i in 1:length(pr)){
      suppressWarnings({
        err = rbind(err, mdta[, .(id = c('D', 'p.value', 'alternative', 'method', 'call'), scale = pr[i], unlist(ks.test(eval(parse(text = pr[i]))[RID == 'FROM'], eval(parse(text = pr[i]))[RID == 'TO']) )), by = variable ])
      })}

    cerr = dcast(err, variable + scale ~id, value.var = 'V3')
    DD = cerr[, sum(as.double(D)), by = variable][, V1]
    pp = cerr[, sum(as.numeric(p.value) > .1), by = variable][, V1]

    idif = pdif - DD
    pdif = DD
    ite = ite + 1

  }

  cat('\n')
  list(dta = copy(mdta[, .(variable, RID, DTM, D1 = ORIG)]))


}
