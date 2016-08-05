library(musica)

data("basin_PT")

scen = basin_PT$sim_scen
ctrl = basin_PT$sim_ctrl
obs = basin_PT$obs_ctrl

rev_difv = function(x, v){
  if (v %in% getOption('additive_variables')) {return(sum(x))} else {return(prod(x))}

}


apxfun = function(y){
  function(v) y
}


abb = month
#dta = rbind(data.table(RID = 'SCEN', scen), data.table(RID = 'CTRL', ctrl), data.table(RID = 'OBS', obs), data.table(RID = 'MODIF', obs))
#dta = rbind(data.table(RID = 'TO', obs), data.table(RID = 'FROM', ctrl), data.table(RID = 'NEWDATA', scen), data.table(RID = 'MODIF', scen))
dta = rbind(data.table(RID = 'TO', obs), data.table(RID = 'FROM', ctrl), data.table(RID = 'NEWDATA', ctrl), data.table(RID = 'MODIF', ctrl))
mdta = melt(dta, id.vars = c('RID', 'DTM'), value.name = 'ORIG')

mdta[, G1 := mean(ORIG), by = .(RID, variable)]
mdta[, Y1 := dif(mean(ORIG), G1, variable[1]), by = .(RID, variable, year(DTM))]
mdta[, M1 := dif(mean(ORIG), rev_dif(G1, Y1, variable[1]), variable[1]), by = .(RID, variable, year(DTM), month(DTM))]
mdta[, D1 := dif(ORIG, rev_difv(c(G1, Y1, M1), variable[1]), variable[1]), by = .(RID, variable, DTM)]
mdta[, GRP := month(DTM)]
#mdta[, C := rev_difv(c(G1, Y1, M1, D1), variable[1]), by = .(RID, variable, DTM)]

pr = c('G1', 'Y1', 'M1', 'D1')
cpr = code2period(pr)
repor = list()
plots = list()
pdif = c(NA, NA)
cdif = c(NA, NA)
idif = pdif - cdif
ite = 1
plots[[ite]] = list()


i = 4
wet_int_thr = .1
ip = seq(0, 1, by = .01)
model = 'identity'
sfce_ = provide_sfce(model)

for (i in 1:length(pr)){

  mdta[, MEA := eval(parse(text = pr[i]))]
  mdta[, p := NA_real_]
  mdta[, cut := as.IDate(musica:::cut.Date(DTM, breaks = cpr[i])), by = .(RID, variable)]
  mdta[, GRP := month(cut)]

  mdta[(variable == 'TAS') | (variable=='PR' & MEA > wet_int_thr), p := prob(MEA), by = .(RID, variable, GRP)]

  smry = mdta[!is.na(p), .(MEA = MEA[1], p = p[1]), by = .(RID, variable, cut, GRP)]

  smry = if (cpr[i]==0){
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

  plots[[ite]][[pr[i]]] = ggplot(cf) + geom_line(aes(x = p, y = CF, group = GRP)) + facet_wrap(GRP~variable, scale = 'free')

  f = list()
  for (g in cf[, unique(GRP)]){
    for (v in cf[, unique(variable)]){
      f[[length(f)+1]] = if (cpr[i]!=0) {
        cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(approxfun(x = p, y = CF, rule = 2)))]
      } else {
        cf[GRP == g & variable == v, .(GRP = g, variable = v, fun = .(apxfun( CF )))]
      }

    }
  }

  f = rbindlist(f)

  mdta = f[mdta, on = c('GRP', 'variable')]
  mdta[, CF:= fun[1][[1]](p), by = .(variable, GRP) ]
  mdta[RID=='MODIF', pr[i] := rev_dif(eval(parse(text = pr[i])), CF, variable[1]), by = .(GRP, variable)]
  #dta[is.na(D1), D1:=0]
  mdta[, c('GRP', 'fun', 'cut', 'MEA', 'p', 'CF') := NULL]

}

mdta[is.na(D1), D1 := 0]
mdta[, ORIG := rev_difv(c(G1, Y1, M1, D1), variable[1]), by = .(RID, variable, DTM)]
