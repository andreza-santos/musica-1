library(bilan)
#library(gr4j)
library(musica)
library(parallel)

scales = rev(c('D1' = '1 day', 'D10' = '10 days', 'M1' = '1 month', 'M3' = '3 months', 'M6' = '6 months', 'Y1' = '1 year', 'Y5' = '5 year'))
ys = -months(1)
ab = quarter


pov = 'CUCICE'
setwd('/home/hanel/LAPV/paper_bc/data/bias_corrected/')
dta = readRDS(paste0(pov, '.rds'))

a = dta[, .(strsplit(SID, '_')), by = SID]
a[, EXP:=sapply(V1, function(x)x[3])]
a[, GCM:=sapply(V1, function(x)x[2])]
a[, RCM:=sapply(V1, function(x)x[5])]
a[, RUN:=sapply(V1, function(x)x[4])]
aa = a[EXP!='historical', .(SID = unique(SID), SUB_EXP = c('historical', EXP)), by = .(GCM, RCM, RUN, EXP)]
aa[, sid := (1:.N+1) %/% 2 ]
aa[, RID := paste(GCM, RCM, EXP, RUN, sep = '_')]
aa[SUB_EXP =='historical', SID:= gsub('rcp26|rcp45|rcp60|rcp85', 'historical', SID)]


setwd('/home/hanel/LAPV/paper_bc/data/bil_kal/')
.b = bil.new(type = 'd', file = paste0(pov, '.bil'))

o = data.table(dcast(dta[year(DTM) %in% c(1970:1999) & ID == 'd' & SID == SID[1], .(DTM, variable, OBS)], DTM ~ variable))

.o = data.table(o, run_bil(o, vars = 'RM', .b))
obs = decomp(.o[, .(DTM, PR, TAS, RM)], scale = scales, year_starts = ys, agg_by = ab)
FU = list(Q90 = Q(.9), Q50 = Q(.5), SD = sd, MEAN = mean)



for (dy in c(3, 5, 10, 30)){

  #for (i in 1:500) {
  gc()
  RR = mclapply(1:500, function(i){
    R = list()
    cat(dy, i, '\n', sep = '\t')

    if (i < 38){
      dt = dta[ SID %in% aa[sid == i, SID] & ID == 'd']
      s = data.table(dcast(dt[year(DTM) %in% c(1970:1999) , .(DTM, variable, SIM)], DTM ~ variable))
      setnames(s, c('PR', 'TAS'), c('sPR', 'sTAS'))

      os = o[s, on = 'DTM']
    } else {os = copy(o)}

    os[, D := cut_width((DTM), dy)]
    os[, mD:= cut_width(yday(DTM), dy)]
    os[, d := sample(1:.N), by = D]
    os[, md:= sample(1:.N), by = mD]

    o1 = os[, .(ID = 'within', DTM, PR = as.double(sample(PR)), TAS = as.double(sample(TAS))), by = D]
    o2 = os[, .(ID = 'within_dep', DTM, PR = PR[d], TAS = TAS[d]), by = D]
    o3 = os[, .(ID = 'out', DTM, PR = as.double(sample(PR)), TAS = as.double(sample(TAS))), by = mD]
    setkey(o3, DTM)
    o4 = os[, .(ID = 'out_dep', DTM, PR = PR[md], TAS = TAS[md]), by = mD]
    setkey(o4, DTM)
    if (i < 38){
      o5 = os[, .(ID = 'assim', DTM, PR = sort(PR)[rank(sPR)],  TAS = sort(TAS)[rank(sTAS)]), by = D]
      o6 = os[, .(ID = 'assim_out', DTM, PR = sort(PR)[rank(sPR)],  TAS = sort(TAS)[rank(sTAS)]), by = mD]
      setkey(o6, DTM)
    }
    O = rbind(o1[, .(ID, DTM, PR, TAS)], o2[, .(ID, DTM, PR, TAS)], o3[, .(ID, DTM, PR, TAS)], o4[, .(ID, DTM, PR, TAS)])
    if (i < 38) O = rbind(O, o5[, .(ID, DTM, PR, TAS)], o6[, .(ID, DTM, PR, TAS)])
    O = data.table(O, RM = O[, run_bil(data.table(DTM, PR, TAS), vars = 'RM', b = .b), by = ID]$RM)
    so = split(O, O$ID)
    do = lapply(so, function(x)decomp(x[,.(DTM, PR, TAS, RM)], scales = scales, year_starts = ys, agg_by = ab))

    for (fu in 1:length(FU)){
      bi = compare(do, compare_to = obs, fun = FU[[fu]])
      R[[length(R)+1]] = data.table(I = i, FUN = names(FU)[fu], bi, SID = if (i<60) {aa[sid == i, SID[1]]} else {NA})

      }

    return(R)
    }, mc.cores = 8)
  print(dy)
  saveRDS(RR, gsub('XXX', dy, '/home/hanel/GACR/BC/resamp_XXX.rds'))
}


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
