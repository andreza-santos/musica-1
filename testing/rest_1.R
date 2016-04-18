library(bilan)
library(gr4j)
library(musica)

pov = 'CUCICE'
setwd('/home/hanel/LAPV/paper_bc/data/bias_corrected/')
dta = readRDS(paste0(pov, '.rds'))
dta[, RID:= gsub('historical|rcp26|rcp45|rcp60|rcp85', '', SID), by = SID]

setwd('/home/hanel/LAPV/paper_bc/data/bil_kal/')
.b = bil.new(type = 'd', file = paste0(pov, '.bil'))

o = data.table(dcast(dta[year(DTM) %in% c(1970:1999) & ID == 'd' & RID == RID[1], .(DTM, variable, OBS)], DTM ~ variable))
s = data.table(dcast(dta[year(DTM) %in% c(1970:1999) & ID == 'd' & RID == RID[1], .(DTM, variable, SIM)], DTM ~ variable))

# kontrolni korigovana -  bude je v datasetu, nebo znovu koriguj
#k = data.table(dcast(dta[year(DTM) %in% c(1970:1999) & ID == 'd' & RID == RID[1], .(DTM, variable, COR)], DTM ~ variable))#trans(o, s, fun = correct)
k = correct(o, s)

# budouci simulace - nutno osetrit, aby tam nespadly ruzne RCP scenare, tj. nrow(f) by melo byt zhruba 10950
f = data.table(dcast(dta[year(DTM) %in% c(2070:2999) & ID == 'd' & RID == RID[1], .(DTM, variable, SIM)], DTM ~ variable))

# budouci korigovana - bude je v datasetu, nebo znovu koriguj
#cf = data.table(dcast(dta[year(DTM) %in% c(2070:2999) & ID == 'd' & RID == RID[1], .(DTM, variable, COR)], DTM ~ variable))
cf = correct(o, s, f)

o = data.table(o, run_bil(o, vars = 'RM',  b = .b))
s = data.table(s, run_bil(s, vars = 'RM',  b = .b))
k = data.table(k, run_bil(k, vars = 'RM',  b = .b))
f = data.table(f, run_bil(f, vars = 'RM',  b = .b))
cf = data.table(cf, run_bil(cf, vars = 'RM',  b = .b))

# definuje dekompozici
scales = rev(c('D1' = '1 day', 'D10' = '10 days', 'M1' = '1 month', 'M3' = '3 months', 'M6' = '6 months', 'Y1' = '1 year', 'Y5' = '5 year'))
# kdy zacina rok?
ys = -months(1)
# agregovat po mesicich (ab = month) nebo sezonach (ab = quarter)?
ab = quarter

# dekompozice
obs = decomp(o, scales = scales, year_starts = ys, agg_by = ab)
sim = decomp(s, scales = scales, year_starts = ys, agg_by = ab)
kor = decomp(k, scales = scales, year_starts = ys, agg_by = ab)
fsim =decomp(f, scales = scales, year_starts = ys, agg_by = ab)
fkor =decomp(cf, scales = scales, year_starts = ys, agg_by = ab)

# porovnani
F = median#Q(.9)
bi = compare(list(RESIDUAL_BIAS = kor), compare_to = obs, fun = F, wet_int_only = TRUE, exclude_below = .9)
bi2 = compare(list(DIF_CHANGE = fkor), compare_to = fsim, fun = F, wet_int_only = TRUE, exclude_below = .9)
bi = rbind(bi, bi2)

ggplot(bi) + geom_line(aes(x = TS, y = DIF, col = sscale2sea(sub_scale, ys), group = sub_scale)) + facet_wrap(variable~comp, scale = 'free', ncol = 2)+scale_x_log10(breaks = tscale(names(scales)), lab = names(scales)) + theme_linedraw() + theme(legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = element_rect(fill = NA))+scale_color_brewer(palette = 'Set1', name ='SEASON')

# korelace

bi = vcompare(list(OBS = obs, SIM = sim, COR = kor), fun = cor)
bi[, SEA:=sscale2sea(sub_scale, ys)]
ggplot(bi) + geom_line(aes(x = TS, y = value, col = ID))+facet_grid(VARS~SEA, scales = 'free') + scale_x_log10(breaks = tscale(names(scales)), lab = names(scales))+ theme_linedraw() + theme(legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = element_rect(fill = NA))+scale_color_brewer(palette = 'Set1', name ='SIMULATION')
