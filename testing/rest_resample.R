library(bilan)
library(gr4j)
library(musica)

scales = rev(c('D1' = '1 day', 'D10' = '10 days', 'M1' = '1 month', 'M3' = '3 months', 'M6' = '6 months', 'Y1' = '1 year', 'Y5' = '5 year'))
ys = -months(1)
ab = quarter


pov = 'CUCICE'
setwd('/home/hanel/LAPV/paper_bc/data/bias_corrected/')
dta = readRDS(paste0(pov, '.rds'))
dta[, RID:= gsub('historical|rcp26|rcp45|rcp60|rcp85', '', SID), by = SID]

setwd('/home/hanel/LAPV/paper_bc/data/bil_kal/')
.b = bil.new(type = 'd', file = paste0(pov, '.bil'))

o = data.table(dcast(dta[year(DTM) %in% c(1970:1999) & ID == 'd' & RID == RID[1], .(DTM, variable, OBS)], DTM ~ variable))
s = data.table(dcast(dta[year(DTM) %in% c(1970:1999) & ID == 'd' & RID == RID[1], .(DTM, variable, SIM)], DTM ~ variable))
setnames(s, c('PR', 'TAS'), c('sPR', 'sTAS'))

os = o[s, on = 'DTM']

o = data.table(o, run_bil(o, vars = 'RM', .b))
obs = decomp(o[, .(DTM, RM)], scale = scales, year_starts = ys, agg_by = ab)

dy = 10

  os[, D := cut(DTM, breaks = paste(dy, 'days'))]
  os[, d := sample(1:.N), by = D]
  os[, md:= sample(1:.N), by = month(D)]

  o1 = os[, .(ID = 'within', DTM, PR = sample(PR), TAS = sample(TAS)), by = D]
  o2 = os[, .(ID = 'within_dep', DTM, PR = PR[d], TAS = TAS[d]), by = D]
  o3 = os[, .(ID = 'out', DTM, PR = sample(PR), TAS = sample(TAS)), by = month(D)]
  setkey(o3, DTM)
  o4 = os[, .(ID = 'out_dep', DTM, PR = PR[md], TAS = TAS[md]), by = month(D)]
  setkey(o4, DTM)
  o5 = os[, .(ID = 'assim', DTM, PR = sort(PR)[rank(sPR)],  TAS = sort(TAS)[rank(sTAS)]), by = D]
  o6 = os[, .(ID = 'assim_out', DTM, PR = sort(PR)[rank(sPR)],  TAS = sort(TAS)[rank(sTAS)]), by = month(D)]
  setkey(o6, DTM)
  O = rbind(o1[, .(ID, DTM, PR, TAS)], o2[, .(ID, DTM, PR, TAS)], o3[, .(ID, DTM, PR, TAS)], o4[, .(ID, DTM, PR, TAS)], o5[, .(ID, DTM, PR, TAS)], o6[, .(ID, DTM, PR, TAS)])
  O = data.table(O, RM = O[, run_bil(data.table(DTM, PR, TAS), vars = 'RM', b = .b), by = ID]$RM)
  so = split(O, O$ID)
  do = lapply(so, function(x)decomp(x[,.(DTM, RM)], scales = scales, year_starts = ys, agg_by = ab))

  bi = compare(do, compare_to = obs, fun = Q(.9))



# porovnani
F = Q(.9)

bi = compare(list(RESIDUAL_BIAS = kor), compare_to = obs, fun = F, wet_int_only = TRUE, exclude_below = .9)
bi2 = compare(list(DIF_CHANGE = fkor), compare_to = fsim, fun = F, wet_int_only = TRUE, exclude_below = .9)
bi = rbind(bi, bi2)

ggplot(bi) + geom_line(aes(x = TS, y = DIF, col = sscale2sea(sub_scale, ys), group = sub_scale)) + facet_wrap(variable~comp, scale = 'free', ncol = 2)+scale_x_log10(breaks = tscale(names(scales)), lab = names(scales)) + theme_linedraw() + theme(legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = element_rect(fill = NA))+scale_color_brewer(palette = 'Set1', name ='SEASON')

# korelace

bi = vcompare(list(OBS = obs, SIM = sim, COR = kor), fun = cor)
bi[, SEA:=sscale2sea(sub_scale, ys)]
ggplot(bi) + geom_line(aes(x = TS, y = value, col = ID))+facet_grid(VARS~SEA, scales = 'free') + scale_x_log10(breaks = tscale(names(scales)), lab = names(scales))+ theme_linedraw() + theme(legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = element_rect(fill = NA))+scale_color_brewer(palette = 'Set1', name ='SIMULATION')
