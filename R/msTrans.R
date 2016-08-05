# #' @export
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
#       return(exp(sfce(log(y), x, v)))
#     }
#
#   }
#   sfce_
# }

#' @export
ms_trans = function(dta, model = 'const', model_par = list(NULL), agg_by = month, wet_int_thr = 0.1, maxiter = 10, tol = 1e-4, report = FALSE, scale = c('G1', 'Y1', 'M1', 'D1'), ipreci = .001, year_starts = months(0)){

  apxfun = function(y){
    function(v) y
  }

  pr = code2period(scale)
  ip = seq(0, 1, by = ipreci)

  repor = list()
  plots = list()
  pdif = c(NA, NA)
  cdif = c(NA, NA)
  idif = pdif - cdif
  ite = 1

  sfce = if (!is.function(model)){
    switch(model,
           'loess' = function(y, x, v){loess(y ~ x, model_par[[v]])$fitted},
           'lm' = function(y, x, v){lm(y ~ x, model_par[[v]])$fitted},
           'const' = function(y, x, v){ rep_len(mean(y), length(y)) },#lm(y ~ 1, model_par[[v]])$fitted},
           'smooth' = function(y, x, v){c(smooth(y, model_par[[v]]))},
           'runmed' = function(y, x, v){runmed(y, model_par[[v]])},
           'smooth.spline' = function(y, x, v){smooth.spline(x, y, model_par[[v]])$y},
           'identity' = function(y, x, v){y}
    )
  } else {function(...)model(x, y, model_par[[v]])}

  sfce_ = function(y, x, v){

    if (length(x)==1) return(y)
    if (v %in% getOption('additive_variables')) {
      return(sfce(y, x, v))
    } else {
      return((sfce((y), x, v))^1)
    }

  }


  while ( (ite <= maxiter) & ((max(abs(idif)) > tol) | any(is.na(idif))) ) {

    cat('\n', format(ite, width =3), '(', format(signif(idif, 5), width = 12), ') :', '\t')
    plots[[ite]] = list()

    i =1
    for (i in 1:length(pr)){

      cat(names(pr[i]), '\t')
      dta[, cut := as.IDate(musica:::cut.Date(DTM, breaks = pr[i])), by = .(RID, variable)]
      dta[, MEA := mean(D1), by = .(RID, variable, cut)]
      dta[, GRP := agg_by(cut %m-% year_starts)]
      dta[, p := NA_real_]
      dta[(variable == 'TAS') | (variable=='PR' & D1 > wet_int_thr), p := prob(MEA), by = .(RID, variable, GRP)]

      smry = dta[!is.na(p), .(MEA = mean(D1), p = p[1]), by = .(RID, variable, cut, GRP)]

      smry = if (pr[i]==0){
        smry[, .(RID, variable, GRP, p, value = MEA)]
      } else {
        smry[, .(p = ip, value = approx(x = p, xout = ip, y = MEA, rule = 2)$y), by = .(RID, variable, GRP)]
      }

      cf = dcast.data.table(smry, variable + GRP + p ~ RID)
      cf[, sFROM := sfce_(y = FROM, x = p, v = as.character(variable[1])), by = .(variable, GRP)]
      cf[, sMODIF := sfce_(y = MODIF, x = p, v = as.character(variable[1])), by = .(variable, GRP)]
      cf[, sNEWDATA := sfce_(y = NEWDATA, x = p, v = as.character(variable[1])), by = .(variable, GRP)]
      cf[, sTO := sfce_(y = TO, x = p, v = as.character(variable[1])), by = .(variable, GRP)]

      cf[, c('REQ', 'CURR') := .(dif(sTO, sFROM, variable[1]), dif(sMODIF, sNEWDATA, variable[1])), by = .(variable)]
      cf[, CF := dif(REQ, CURR, variable[1]), by = variable]

      plots[[ite]][[names(pr[i]) ]] = ggplot(cf) + geom_line(aes(x = p, y = CF, group = GRP)) + facet_wrap(GRP~variable, scale = 'free')

      f = list()
      for (g in cf[, unique(GRP)]){
        for (v in cf[, unique(variable)]){
          f[[length(f)+1]] = if (pr[i]!=0) {
            cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(approxfun(x = p, y = CF, rule = 2)))]
          } else {
            cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(apxfun( CF )))]
          }

        }
      }

      f = do.call(rbind, f)

      dta = f[dta, on = c('GRP', 'variable')]
      dta[, CF:= fun[1][[1]](p), by = .(variable, GRP) ]
      dta[RID=='MODIF', D1 := rev_dif(D1, CF, variable[1]), by = .(GRP, variable)]
      dta[is.na(D1), D1:=0]
      dta[, c('GRP', 'fun', 'cut', 'MEA', 'p', 'CF') := NULL]

      repor[[length(repor)+1]] = data.table(scale = names(pr)[i], iter = ite, cf)

    }

    r = do.call(rbind, repor)
    cdif = r[iter==ite, ifelse(variable=='TAS', mean(abs(CF)), mean(abs(CF-1))), by = variable][, V1]
    idif = pdif - cdif
    pdif = cdif
    ite = ite + 1
    if (!report) repor = list()

  }


  pplot = function(iter, scale){
    plots[[iter]][[scale]]
  }

  if (report) r = do.call(rbind, repor)

  cat('\n')
  list(dta = copy(dta), pplot = pplot, report = if (report) {r} else {NULL})

}

ext_cast = function(out, what){
  dcast.data.table(out$dta[RID==what, .(variable, DTM, D1)], DTM ~ variable, value.var = 'D1')
}

#
# ################################################
# #' @export
# ms_atrans = function(dta, model = 'const', model_par = list(NULL), agg_by = month, wet_int_thr = 0.1, maxiter = 10, tol = 1e-4, report = FALSE, scale = c('G1', 'Y1', 'M3', 'M1', 'D1'), ipreci = .001, year_starts = months(0)){
#
#   apxfun = function(y){
#     function(v) y
#   }
#
#   deco = dta[, decomp(data.table(DTM, PR, TAS), periods = code2period(scale)[scale!='G1'], full_return = TRUE, agg_by = month), by = RID]
#   cdeco = dcast.data.table(deco[, .(RID, variable, DTM = dDTM, period, value, sub_period)], RID+variable+DTM+sub_period ~ period)
#
#   pe = code2period(scale)
#   ip = seq(0, 1, by = ipreci)
#
#   repor = list()
#   plots = list()
#   pdif = c(NA, NA)
#   cdif = c(NA, NA)
#   idif = pdif - cdif
#   ite = 1
#
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
#       return(exp(sfce(log(y), x, v)))
#     }
#
#   }
#
#
#   while ( (ite <= maxiter) & ((max(abs(idif)) > tol) | any(is.na(idif))) ) {
#
#     cat('\n', format(ite, width =3), '(', format(signif(idif, 5), width = 12), ') :', '\t')
#     plots[[ite]] = list()
#
#     i =1
#
#     for (i in 1:length(pe)){
#
#       cdeco[, cut := as.IDate(musica:::cut.Date(DTM, breaks = pe[i])), by = .(RID, variable)]
#       cdeco[, sub_period := agg_by(cut %m-% year_starts)]
#       cdeco[, MEA := mean(eval(parse(text = names(pe)[max(i, 1)]))), by = .(variable, sub_period, RID, cut)]
#
#       cdeco[, cut := as.IDate(musica:::cut.Date(DTM, breaks = pe[max(i-1, 1)])), by = .(RID, variable)]
#
#       #if (i>1) {
#         #cdeco[, MEA := mean(eval(parse(text = names(pe)[max(i, 1)]))), by = .(variable, sub_period, RID)]
#         cdeco[, prev_MEA := mean(eval(parse(text = names(pe)[max(i-1, 1)]))), by = .(variable, sub_period, RID, cut)]
#       #} else {#cdeco[, MEA := mean(D1), by = .(variable, sub_period, RID)]
#       #  cdeco[, MEA := mean(eval(parse(text = names(pe)[i]))), by = .(variable, sub_period, RID)]
#       #}
#
#       if (i == 1) {cdeco[, ANOM := G1, by = .(variable, sub_period, RID) ]
#       } else {
#         cdeco[, ANOM:= dif(eval(parse(text = names(pe)[i])), prev_MEA, variable[1]), by = .(variable, sub_period, RID)]
#       }
#       cdeco[, p := NA_real_]
#       cdeco[(variable == 'TAS') | (variable=='PR' & eval(parse(text = names(pe)[i])) > wet_int_thr), p := prob(ANOM), by = .(RID, variable, sub_period)]
#
#
#       #smry = cdeco[!is.na(p), .(MEA = mean(D1), p = p[1]), by = .(RID, variable, cut, GRP)]
#       smry = if (pe[i]==0){
#         cdeco[!is.na(p), .(RID = RID[1], variable = variable[1], sub_period = sub_period[1], p = p[1], value = ANOM[1]), by = .(RID, variable)]
#       } else {
#         cdeco[, .(p = ip, value = approx(x = p, xout = ip, y = ANOM, rule = 2)$y), by = .(RID, variable, sub_period)]
#       }
#
#       cf = dcast.data.table(smry[, .(variable, sub_period, p, RID, value)], variable + sub_period + p ~ RID)
#
#       cf[, sFROM := sfce_(y = FROM, x = p, v = as.character(variable[1])), by = .(variable, sub_period)]
#       cf[, sMODIF := sfce_(y = MODIF, x = p, v = as.character(variable[1])), by = .(variable, sub_period)]
#       cf[, sNEWDATA := sfce_(y = NEWDATA, x = p, v = as.character(variable[1])), by = .(variable, sub_period)]
#       cf[, sTO := sfce_(y = TO, x = p, v = as.character(variable[1])), by = .(variable, sub_period)]
#
#       cf[, c('REQ', 'CURR') := .(dif(sTO, sFROM, variable[1]), dif(sMODIF, sNEWDATA, variable[1])), by = .(variable)]
#       cf[, CF := dif(REQ, CURR, variable[1]), by = variable]
#
#       plots[[ite]][[names(pe[i]) ]] = ggplot(cf) + geom_line(aes(x = p, y = CF, group = sub_period)) + facet_wrap(sub_period~variable, scale = 'free')
#
#       f = list()
#       for (g in cf[, unique(sub_period)]){
#         for (v in cf[, unique(variable)]){
#           f[[length(f)+1]] = if (pe[i]!=0) {
#             cf[sub_period == g & variable == v, .(sub_period = g, variable = v, fun = .(approxfun(x = p, y = CF, rule = 2)))]
#           } else {
#             cf[sub_period == g & variable == v, .(sub_period = g, variable = v, fun = .(apxfun( CF )))]
#           }
#
#         }
#       }
#
#       f = do.call(rbind, f)
#
#       cdeco = f[cdeco, on = c('sub_period', 'variable')]
#       cdeco[, CF:= fun[1][[1]](p), by = .(variable, sub_period) ]
#       cdeco[RID=='MODIF', cANOM := rev_dif(ANOM, CF, variable[1]), by = .(sub_period, variable)]
#       #cdeco[RID=='MODIF', names(pe[i]) := cANOM]#rev_dif(cANOM, MEA, variable[1]), by = .(sub_period, variable)]
#       if (i!=1) {
#         cdeco[RID=='MODIF', names(pe[i]) := rev_dif(cANOM, prev_MEA, variable[1]), by = .(sub_period, variable)]
#       } else {
#         cdeco[RID=='MODIF', names(pe[i]) := cANOM, by = .(sub_period, variable)]
#         }
#
#       cdeco[is.na(eval(parse(text = names(pe[i]) ))), names(pe[i]) := 0]
#       #cdeco[is.na(D1), D1 := 0]
#
#       cdeco[, c('fun', 'cut', 'MEA', 'p', 'CF', 'ANOM', 'cANOM') := NULL]
#
#       repor[[length(repor)+1]] = data.table(scale = names(pe)[i], iter = ite, cf)
#
#     }
#
#     # mdeco = melt(cdeco, id.vars = c('sub_period', 'variable', 'RID', 'DTM'), value.name = 'D1', variable.name = 'period')
#     # res  = mdeco[RID=='MODIF', ifelse(unique(variable) == 'PR', prod(D1), sum(D1)), by = .(sub_period, variable, RID, DTM)]
#     # cdeco[RID=='MODIF', D1 := res$V1]
#     #
#     # rdeco = cdeco[RID=='MODIF', decomp(dcast.data.table(data.table(DTM, variable, D1), DTM ~ variable, value.var = 'D1'), periods = code2period(scale)[scale!='G1'], full_return = TRUE, agg_by = month), by = RID]
#     # rcdeco = dcast.data.table(rdeco[, .(RID, variable, DTM = dDTM, period, value, sub_period)], RID+variable+DTM+sub_period ~ period)
#     #
#     # for (ii in 1:length(scale)){
#     #   cdeco[RID=='MODIF', scale[ii] := rcdeco[, .(eval(parse(text = scale[ii]))) ] ]
#     # }
#
#     r = do.call(rbind, repor)
#     cdif = r[iter==ite, ifelse(variable=='TAS', mean(abs(CF)), mean(abs(CF-1))), by = variable][, V1]
#     idif = pdif - cdif
#     pdif = cdif
#     ite = ite + 1
#     if (!report) repor = list()
#
#   }
#
#
#   pplot = function(iter, scale){
#     plots[[iter]][[scale]]
#   }
#
#   if (report) r = do.call(rbind, repor)
#
#   cat('\n')
#   list(dta = copy(cdeco), pplot = pplot, report = if (report) {r} else {NULL})
#
# }
