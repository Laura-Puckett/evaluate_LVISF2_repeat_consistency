library(ggplot2)
inDir = '//nau.froot.nau.edu/cirrus/scratch/lkp58/LVIS/accuracy_assessment/output/'
# files = list.files(inDir, pattern = '.csv')
meta = read.csv('//nau.froot.nau.edu/cirrus/scratch/lkp58/LVIS/accuracy_assessment/inputs.csv') %>%
  dplyr::select(outname, Site, TimeDiff, Burned, Status)
meta = as.data.frame(meta, stringsAsFactors=F)

meta = meta[meta$Status %in% c("COMPLETED","ALMOST"),]

combined = data.frame()
for(i in 1:nrow(meta)){
  print(meta[i,"outname"])
  temp =  read.csv(paste0(inDir, meta[i,"outname"], '.csv'))
  
  if("SENSITIVITY_2" %in% colnames(temp)==F){
    temp = cbind("SENSITIVITY_2" = NA, temp)
  }
  if("SENSITIVITY_1" %in% colnames(temp)==F){
    temp = cbind("SENSITIVITY_1" = NA, temp)
  }
  
  temp = temp %>% select(ZG_diff, 
                         SENSITIVITY_1,
                         SENSITIVITY_2,
                         AZIMUTH_1,
                         AZIMUTH_2,
                         INCIDENTANGLE_1,
                         INCIDENTANGLE_2,
                         GLAT_1,
                         CHAN_ZG,
                         RH50_1,
                         RH98_1,
                         RH98_diff,
                         distance)
  
  temp = cbind("name" = meta[i, "outname"],
               "Site" = meta[i, "Site"],
               "TimeDiff" = meta[i, "TimeDiff"],
               "Burned" = meta[i, "Burned"],
               temp)
  combined = rbind(combined, temp)
  
  # print(quantile(temp$ZG_diff, prob = c(0.025, .975), na.rm=T))
  # print(quantile(temp$ZG_diff, prob = c(0.05, .95), na.rm=T))
  # print(quantile(temp$ZG_diff, prob = c(0.25, 0.50, .75), na.rm=T))
}

combined = as.data.frame(combined)
combined = combined %>% dplyr::mutate(SENS_min = ifelse(as.numeric(SENSITIVITY_1)> as.numeric(SENSITIVITY_2),SENSITIVITY_2,SENSITIVITY_1),
                                      AZIM_diff = AZIMUTH_2 - AZIMUTH_1,
                                      INCI_max = ifelse(INCIDENTANGLE_1 >INCIDENTANGLE_2, INCIDENTANGLE_1, INCIDENTANGLE_2))

# pdf(paste0('//nau.froot.nau.edu/cirrus/scratch/lkp58/LVIS/accuracy_assessment/output/',trial,'.pdf'))

p = ggplot(combined, aes(x=ZG_diff)) + 
  geom_histogram(breaks = seq(-1,1,.01), color="black", fill="white")
print(p)

print("same day subset")
subset_sameDay = combined %>% filter(as.character(TimeDiff) %in% c("same"))
print(quantile(subset_sameDay$ZG_diff, prob = c(0.025, 0.25, 0.5, 0.75, .975), na.rm=T))
p = ggplot(subset_sameDay, aes(x=ZG_diff)) + 
  geom_histogram(breaks = seq(-1,1,.01), color="black", aes(fill = Site))
print(p)

print("diff Year")
subset_diffYear = combined %>% filter(as.character(TimeDiff) %in% c("years"))
print(quantile(subset_diffYear$ZG_diff, prob = c(0.025, 0.25, 0.5, 0.75, .975), na.rm=T))
p = ggplot(subset_diffYear, aes(x=ZG_diff)) + 
  geom_histogram(breaks = seq(-1,1,.01), color="black", aes(fill = as.character(Burned)))
print(p)

print("burned subset")
subset_burned = combined %>% filter(as.character(Burned) == "yes")
print(quantile(subset_burned$ZG_diff, prob = c(0.025, 0.25, 0.5, 0.75, .975), na.rm=T))
p = ggplot(subset_burned, aes(x=ZG_diff)) + 
  geom_histogram(breaks = seq(-1,1,.01), color="black", aes(fill = as.character(Site)))
print(p)


print(quantile(combined$ZG_diff, prob = c(0.025, .975), na.rm=T))
print(quantile(combined$ZG_diff, prob = c(0.05, .95), na.rm=T))
print(quantile(combined$ZG_diff, prob = c(0.25, 0.50, .75), na.rm=T))

corrgram(combined %>% dplyr::select(ZG_diff,
                                    SENS_min,
                                    AZIM_diff,
                                    INCI_max,
                                    CHAN_ZG,
                                    distance,
                                    RH50_1,
                                    RH98_1,
                                    RH98_diff,
                                    Site,
                                    TimeDiff,
                                    Burned),
         order=NULL,
         lower.panel=panel.shade,
         text.panel=panel.txt,
         upper.panel=panel.pts,
         main="Correlogram")

dev.off()