dif = function(x, y, variable){
  if (variable == 'TAS') return(x-y)
  r = x/y
  r[y==0] = 0
  return(r)
}

rev_dif = function(x, y, variable){
  if (variable == 'TAS') return(x+y)
  r = x*y
  #r[y==0] = 0
  return(r)
}


decomp = function(x, year_starts = months(0), scales = c(Y1 = '1 year', M6 = '6 month', M3 = '3 months', M1 = '1 months', D5 = '15 days', D1 = '1 day'), agg_by = quarter){

  x = copy(x)
  # shift for arbitrary year start (climatological, hydrological etc.)
  x[, cDTM := as.IDate(DTM %m-% year_starts)]

  for (i in 1:length(scales)){
    bby = if (i>1) {names(scales)[1:(i-1)]} else {NULL}
    x[, (names(scales)[i]) := cut(cDTM, breaks = scales[i]), by = bby]#eval(parse(text = bby))]
  }

  wh = names(x)[! names(x) %in% c('DTM', 'cDTM', names(scales))]
  x[, G1:= mean(cDTM)]
  mo = melt(x, measure.vars = wh)

  R = list()
  for (i in 0:length(scales)){
    bby = if (i>0) {names(scales)[i]} else {'G1'}
    R[[i+1]] = mo[,  .(scale = bby, value = mean(value), sub_scale = unique(agg_by(cDTM)), scale_pos = unique(prse), N = .N), by = .(variable, prse = eval(parse(text = bby)))]

    }

  R = do.call(rbind, R)
  R[, TS:= tscale(scale)]
  setnames(R, 'prse', 'DTM')
  R[, DTM := as.IDate(DTM %m+% year_starts)]
  R[, scale_pos := as.IDate(scale_pos %m+% year_starts)]
  #R[, excl:= N < 0.9 * (TS/24)]
  copy(R)

}

tscale = function(x, nyears = 30){

  num = suppressWarnings(as.integer(gsub("[^\\d]+", "", x, perl = TRUE)))
  xx = gsub('G', 'year ',  x)
  xx = gsub('Y', 'year ', xx)
  xx = gsub('M', 'month ', xx)
  xx = gsub('D', 'day ', xx)
  xx[grepl('G', x)] = gsub(1, nyears, xx[grepl('G', x)])

  xs = strsplit(xx, ' ')
  xs = data.table(x, do.call(rbind, xs))
  xs[, ts:= unclass(duration(as.double(V2[1]), V1[1]))/3600, by = x]
  return(xs$ts)

}

compare = function(x, compare_to, fun = mean, wet_int_only = TRUE, wet_int_thr = 0.1, exclude_below = 0.9){

  lst = c(x, COMPARE_TO = list(compare_to))
  for (i in 1:length(lst)){
    lst[[i]][, ID:=names(lst)[i]]
  }
  lst = do.call(rbind, lst)
  lst = lst[N >= exclude_below * (TS/24)]
  stat =
    if (wet_int_only){
      lst[!(variable!='TAS' & value<= wet_int_thr), .(value = fun(value)), by = .(variable, scale, TS, sub_scale, ID)]
    } else {
      lst[, .(value = fun(value)), by = .(variable, scale, TS, sub_scale, ID)]
    }
  cstat = dcast(stat, variable + scale + TS + sub_scale ~ ID)
  stat = melt(cstat, measure.vars = names(x), variable.name = 'comp')
  stat[, .(DIF = dif(value, COMPARE_TO, variable[1])), by = .(variable, scale, TS, sub_scale, comp)]

}

vcompare = function(x, fun = cor, wet_int_only = TRUE, wet_int_thr = 0.1, exclude_below = 0.9){

  lst = x #c(x, COMPARE_TO = list(compare_to))
  for (i in 1:length(lst)){
    lst[[i]][, ID:=names(lst)[i]]
  }
  lst = do.call(rbind, lst)

  vars = lst[, unique(variable)]
  cmb = t(combn(vars, 2))
  map = data.table(TYP = rep(names(x),  each = nrow(cmb)), cmb)

  B = list()
  for (i in 1:nrow(map)){
    yy = dcast(lst[ID == map[i, TYP] & N >= exclude_below * (TS/24), ], ... ~ variable, value.var = 'value')
    cy = yy[, .(ID = map[i, TYP], V1 = map[i, V1], V2 = map[i, V2], TS = tscale(scale), value = fun(eval(parse(text = map[i, V1])),  eval(parse(text = map[i, V2]))), VARS = map[i, paste(V1, V2, sep = ' x ')]), by = .(scale, sub_scale) ]
    B[[i]] = cy
  }

  do.call(rbind, B)

}



Q = function(p, ...){
  function(x)quantile(x, p, ...)
}


run_bil = function(dta, b = NULL, vars = 'RM'){
  if (is.null(b)) {b = bil.new(type = 'd', data = dta[, .(DTM, P = PR, T = TAS)]) } else {
    bil.set.values(b, dta[, .(DTM, P = PR, T = TAS)])
  }
  bil.pet(b)
  d = data.frame(data.frame(bil.run(b))[, vars])
  names(d) = vars
  d
}

correct = function(x, y, z = copy(y), fun = correct, ...){

  mx = melt(x, id.vars = 'DTM', value.name = 'obs')
  my = melt(y, id.vars = 'DTM', value.name = 'sim')
  mz = melt(z, id.vars = 'DTM', value.name = 'pred')
  m = mx[my, on = c('DTM', 'variable')]

  mon = 1
  v = 'PR'
  for (v in m[, unique(variable)]){
    for (mon in 1:12){
      f = m[variable==v & month(DTM) == mon, fitQmap(obs, sim, method = 'QUANT', wet.day = v == 'PR', qstep = .001)]
      mz[variable==v & month(DTM) == mon , cor:=doQmap(pred, f)]
    }
  }

  copy(dcast(mz[, .(DTM, variable, cor)], DTM ~ variable, value.var = 'cor'))
}

# correct = function(x, y, z, var, ...){
#   fi = fitQmap(x, y, method= 'QUANT', wet.day = var == 'PR', ...)
#   doQmap(z, fi, ...)
# }

month2sea = function(dtm, year_starts){
  i = ((1:12-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  i[i==0] = 12
  id = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N' , 'D')[i]
  id = rep(c(paste0(id[1:3], collapse = ''), paste0(id[1:3+3], collapse = ''), paste0(id[1:3+6], collapse = ''), paste0(id[1:3+9], collapse = '')), each = 3)
  i = month(dtm %m-% year_starts)#((month(dtm)-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  i[i==0] = 12
  id[i]
}

sscale2sea = function(sub_scale, year_starts){
  i = ((1:12-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  i[i==0] = 12
  id = c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N' , 'D')[i]
  id = c(paste0(id[1:3], collapse = ''), paste0(id[1:3+3], collapse = ''), paste0(id[1:3+6], collapse = ''), paste0(id[1:3+9], collapse = ''))
  # i = month(dtm %m-% year_starts)#((month(dtm)-1) + month(as.Date('1970-01-01') + year_starts)) %% 12
  # i[i==0] = 12
  id[sub_scale]
}
