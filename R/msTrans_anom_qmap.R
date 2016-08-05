# provide_sfce = function(model, model_par = list(NULL)){
#   sfce = if (!is.function(model)){
#     switch(model,
#            'loess' = function(y, x, v){loess(y ~ x, model_par[[v]])$fitted},
#            'lm' = function(y, x, v){lm(y ~ x, model_par[[v]])$fitted},
#            'const' = function(y, x, v){lm(y ~ 1, model_par[[v]])$fitted},
#            'smooth' = function(y, x, v){c(smooth(y, model_par[[v]]))},
#            'runmed' = function(y, x, v){runmed(y, model_par[[v]])},
#            'smooth.spline' = function(y, x, v){smooth.spline(x, y, model_par[[v]])$y},
#            'identity' = function(y, x, v){y}
#     )
#   } else {function(...)model(x, y, model_par[[v]])}
#
#   sfce_ = function(y, x, v){
#
#     if (length(x)==1) return(y)
#     if (v %in% getOption('additive_variables')) {
#       return(sfce(y, x, v))
#     } else {
#       return((sfce((y), x, v)))
#     }
#
#   }
#   sfce_
# }

approx = function(x, y = NULL, xout, method = "linear", n = 50, yleft, yright, rule = 1, f = 0, ties = mean){
   # print(unique(x))
   # print(unique(y))
   # print(length(xout))

  if (length(unique(x)) == 1)  {
    return(list(x = xout, y = rep(mean(y), length(xout))))
  } else {
    return(
      stats::approx(x, y, xout, method , n , yleft, yright, rule, f, ties)
    )
  }
}

msTrans = function(dta, agg_by = month, wet_int_thr = 0.1, maxiter = 10, tol = 1e-4, report = FALSE, ipreci = .001, year_starts = months(0)){

  scales = c('G1', 'Y1', 'M1', 'D1')
  apxfun = function(y){
    function(v) y
  }

  cpr = code2period(scales)
  ip = seq(0, 1, by = ipreci)

  repor = list()
  plots = list()
  pdif = c(NA, NA)
  cdif = c(NA, NA)
  pp = c(NA, NA)
  idif = pdif - cdif
  ite = 1

  mdta = melt(dta, id.vars = c('RID', 'DTM'), value.name = 'ORIG')

  while ( (ite <= maxiter) & ((max(abs(idif)) > tol) | any(is.na(idif))) ) {

  cat('\n', format(ite, width =3), '(', format(signif(idif, 5), width = 12), format(pp, width = 6), ') :', '\t')
    plots[[ite]] = list()

    i =1
    mdta[, G1 := mean(ORIG), by = .(RID, variable)]
    mdta[, Y1 := dif(mean(ORIG), G1, variable[1]), by = .(RID, variable, year(DTM))]
    mdta[, M1 := dif(mean(ORIG), rev_dif(G1, Y1, variable[1]), variable[1]), by = .(RID, variable, year(DTM), month(DTM))]
    mdta[, D1 := dif(ORIG, rev_difv(c(G1, Y1, M1), variable[1]), variable[1]), by = .(RID, variable, DTM)]

    for (i in 1:length(pr)){

      switch(
        as.character(i),
             '1' = {mdta[, G1 := mean(ORIG), by = .(RID, variable)]},
             '2' = {mdta[, Y1 := dif(mean(ORIG), G1, variable[1]), by = .(RID, variable, year(DTM))]},
             '3' = {mdta[, M1 := dif(mean(ORIG), rev_dif(G1, Y1, variable[1]), variable[1]), by = .(RID, variable, year(DTM), month(DTM))]},
             '4' = {mdta[, D1 := dif(ORIG, rev_difv(c(G1, Y1, M1), variable[1]), variable[1]), by = .(RID, variable, DTM)]}
             )

      #i =1
      # mdta[, G1 := mean(ORIG), by = .(RID, variable)]
      # mdta[, Y1 := dif(mean(ORIG), G1, variable[1]), by = .(RID, variable, year(DTM))]
      # mdta[, M1 := dif(mean(ORIG), rev_dif(G1, Y1, variable[1]), variable[1]), by = .(RID, variable, year(DTM), month(DTM))]
      # mdta[, D1 := dif(ORIG, rev_difv(c(G1, Y1, M1), variable[1]), variable[1]), by = .(RID, variable, DTM)]

      mdta[, MEA := eval(parse(text = pr[i]))]
     # mdta[, p := NA_real_]
      mdta[, cut := as.IDate(musica:::cut.Date(DTM, breaks = cpr[i])), by = .(RID, variable)]
      mdta[, GRP := month(cut %m-% year_starts)]

      #mdta[(variable == 'TAS') | (variable=='PR' & MEA > wet_int_thr), p := prob(MEA), by = .(RID, variable, GRP)]

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

      # mdta[, CF:= fun[1][[1]](p), by = .(variable, GRP) ]
      # mdta[RID=='MODIF', pr[i] := rev_dif(eval(parse(text = pr[i])), CF, variable[1]), by = .(GRP, variable)]
      mdta[, c('GRP', 'fun', 'cut', 'MEA') := NULL]
      #repor[[length(repor)+1]] = data.table(scale = pr[i], iter = ite, cf)

      #if (i == 4) mdta[is.na(D1), D1 := 0]
      mdta[, ORIG := rev_difv(c(G1, Y1, M1, D1), variable[1]), by = .(RID, variable, DTM)]
      #
      # #
      # mdta[, mean(ORIG), by = .(variable, RID)]
      # mdta[, mean(ORIG), by = .(variable, RID, year(DTM))][, sd(V1), by = .(variable, RID)]
      #
      # mo = mdta[, mean(ORIG), by = .(variable, RID, year(DTM), month(DTM))]
      # ggplot(mo[RID %in% c('FROM', 'TO'), .(V1 = sort(V1), year = 1:length(V1)), by = .(month, variable, RID)]) +  geom_line(aes(x = year, y = V1, col = RID))+  geom_point(aes(x = year, y = V1, col = RID)) + facet_grid(variable ~ month, scale = 'free')
      # #a = mdta[, .(ks.test(G1[RID == 'FROM'], G1[RID == 'TO'])$p.value, ) ]


      }

    mdta[, G1 := mean(ORIG), by = .(RID, variable)]
    mdta[, Y1 := dif(mean(ORIG), G1, variable[1]), by = .(RID, variable, year(DTM))]
    mdta[, M1 := dif(mean(ORIG), rev_dif(G1, Y1, variable[1]), variable[1]), by = .(RID, variable, year(DTM), month(DTM))]
    mdta[, D1 := dif(ORIG, rev_difv(c(G1, Y1, M1), variable[1]), variable[1]), by = .(RID, variable, DTM)]

    err = data.table()
    for (i in 1:length(pr)){
      suppressWarnings({
      err = rbind(err, mdta[, .(id = c('D', 'p.value', 'alternative', 'method', 'call'), scale = pr[i], unlist(ks.test(eval(parse(text = pr[i]))[RID == 'FROM'], eval(parse(text = pr[i]))[RID == 'TO']) )), by = variable ])
    })}

    cerr = dcast(err, variable + scale ~id, value.var = 'V3')
    DD = cerr[, sum(as.double(D)), by = variable][, V1]
    pp = cerr[, sum(as.numeric(p.value) > .1), by = variable][, V1]

    #r = do.call(rbind, repor)
    #cdif = r[iter==ite, ifelse(variable=='TAS', mean(abs(CF)), mean(abs(CF-1))), by = variable][, V1]
    idif = pdif - DD
    pdif = DD
    ite = ite + 1
    #if (!report) repor = list()

    }

  pplot = function(iter, scale){
    plots[[iter]][[scale]]
  }

  if (report) r = do.call(rbind, repor)

  cat('\n')
  list(dta = copy(mdta[, .(variable, RID, DTM, D1 = ORIG)]), pplot = pplot, report = if (report) {r} else {NULL})


}

