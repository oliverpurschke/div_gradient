library(maps)
library(ecoretriever)
#library(mobr)
source('./mobr/R/mobr.R')
source('./mobr/R/mobr_boxplots.R')
source('./scripts/div_functions.R')

dir.create('./data/filtered_data')
dir.create('./results')

# data import ---------------------------------------------------------------

## Purpose: to read in the .txt files of the BCI censues and output a cleaned
## csv file that can then be used to perform calculations on. The filtering 
## performed in this script is for a biodiversity analysis and may not be
## appropriate for analyses targeted at other topics
## Metadata: ~/datasets/CTFSplots/BCI/bci50ha.doc 

print('Filtering and cleaning empirical data, ...')

# bci ---------------------------------------------------------------------

dat = read.table('./data/bci_census7.txt', sep='\t', header=TRUE)

good_data = dat$Status == 'alive' & !is.na(dat$DBH) & 
           !is.na(dat$gx) & !is.na(dat$gy) & 
           dat$Latin != 'Unidentified species' &
           dat$Stem == 'main' 
          
dat = dat[good_data, c('Latin', 'gx', 'gy')]

range(dat$gx)
range(dat$gy)

# filter down to smaller spatial scale
# use central plot of 640 to 400
xmin = (1000 - 640)/2
xmax = 1000 - xmin
ymin = (500 - 400) / 2
ymax = 500 - ymin

dat = dat[dat$gx > xmin & dat$gx < xmax &
          dat$gy > ymin & dat$gy < ymax, ]

write.csv(dat, file='./data/filtered_data/bci_census7_filtered.csv',
          row.names=F)


# scbi --------------------------------------------------------------------

fileprefix = 'SCBI_initial_woody_stem_census_2012'
dat = read.csv(paste('./data/', fileprefix, '.csv', sep=''))
# lump all Crataegus
dat$Latin[dat$Latin == 'Crataegus pruinosa'] = 'Crataegus sp'
# 
good_data = dat$Status == 'alive' & !is.na(dat$DBH) & 
            !is.na(dat$gx) & !is.na(dat$gy) & 
            dat$Latin != 'Acer sp' &             # drop 1 indiv
            dat$Latin != 'Carya sp' &            # drop 80 indiv
            dat$Latin != 'Fraxinus sp' &         # drop 6 indiv
            dat$Latin != 'Quercus sp' &          # drop 7 indiv
            dat$Latin != 'Ulmus sp' &            # drop 95 indiv
            dat$Latin != 'Unidentifed unknown' & # drop 692 indiv
            dat$Stem == 'main'
dat = dat[good_data, c('Latin', 'gx', 'gy')]


write.csv(dat, file=paste('./data/filtered_data/', fileprefix, '_filtered.csv', sep=''), 
          row.names=F)

# huss --------------------------------------------------------------------

fileprefix = '00_HUSS_TREE_2013'
dat = read.csv(paste('./data/', fileprefix, '.csv', sep=''), sep=';', stringsAsFactors=FALSE)

# fix typo
dat$SPECIES[which(dat$SPECIES %in% c("Bu"))] <- "BU"

# give latin names to short names:
dat$SPECIES[which(dat$SPECIES %in% c("AS"))] <- "Populus tremula"
dat$SPECIES[which(dat$SPECIES %in% c("BAH"))] <- "Acer pseudoplatanus"
dat$SPECIES[which(dat$SPECIES %in% c("BU"))] <- "Fagus sylvatica"
dat$SPECIES[which(dat$SPECIES %in% c("BUL"))] <- "Ulmus glabra"
dat$SPECIES[which(dat$SPECIES %in% c("EI"))] <- "Quercus sp." # likely "Quercus petraea"
dat$SPECIES[which(dat$SPECIES %in% c("ES"))] <- "Fraxinus excelsior"
dat$SPECIES[which(dat$SPECIES %in% c("FAH"))] <- "Acer campestre"
dat$SPECIES[which(dat$SPECIES %in% c("HBU"))] <- "Carpinus betulus"
dat$SPECIES[which(dat$SPECIES %in% c("HKB"))] <- "Lonicera sp."
dat$SPECIES[which(dat$SPECIES %in% c("KB"))] <- "Prunus avium"
dat$SPECIES[which(dat$SPECIES %in% c("LI"))] <- "Tilia sp." # likely "Tilia cordata"
dat$SPECIES[which(dat$SPECIES %in% c("PF"))] <- "Euonymus europaeus"
dat$SPECIES[which(dat$SPECIES %in% c("SAH"))] <- "Acer platanoides"
dat$SPECIES[which(dat$SPECIES %in% c("WD"))] <- "Crataegus sp."

# select potentially problematic entries: 
index.PROBLEM <- which(dat$PROBLEM %in% c("1"))
# select the ones that are not really a problem:
index.descPROBLEM1 <- grep("Gegenwinkel*", dat$descPROBLEM)
index.descPROBLEM2 <- grep("BHD*", dat$descPROBLEM)
index.descPROBLEM <- c(index.descPROBLEM1, index.descPROBLEM2)
index.problem2 <- index.PROBLEM[(index.PROBLEM %in% index.descPROBLEM)==FALSE]

dat <- dat[-index.problem2, c('SPECIES', 'EASTING_m', 'NORTHING_m')]
names(dat) <- c('Latin', 'gx', 'gy')

range(dat$gx)
range(dat$gy)

# the huss plot is irregularily shaped -> find maximum rectangular extent

# set new x and y columns for huss data:
#plot(dat$gx, dat$gy)
#identify(dat$gx, dat$gy, plot=TRUE)

#xmin: 1550
#xmax: 11965
#ymin: 3580
#ymax: 1550
#dat[c(1550, 11965, 3580),]

4391445-4390995 # gx
5661862-5661438 # gy

# maximum rectangular extent: 450 x 424 -> cut it to 440 x 400, to fit in 40m grid cells (see scbi & bci analysis):
# To do:cut scbi and bci data to the same extent

xmin = 4390995+5
xmax = 4391445-5
ymin = 5661438+12
ymax = 5661862-12
    
dat = dat[dat$gx > xmin & dat$gx < xmax &
          dat$gy > ymin & dat$gy < ymax, ]

# plot(dat$gx, dat$gy)

write.csv(dat.good, file=paste('./data/filtered_data/', fileprefix, '_filtered.csv', sep=''), row.names=F)


# korup -------------------------------------------------------------------



# gentry ------------------------------------------------------------------

dat = fetch('Gentry')

# examine counts table 
cts = dat$counts

good_data = !is.na(cts$line) &
            !is.na(cts$count) &
            cts$line != 0 & 
            cts$line != 11

cts = cts[good_data, ]

#for only sites with 10 lines
line_cts = tapply(cts$line, list(as.character(cts$site_code)), function(x) length(unique(x)))
gd_sites = names(line_cts[line_cts == 10])
cts = cts[cts$site_code %in% gd_sites, ]

# filter site table down to only new world samples
gd_contintents = c('North America', 'Mesoamerica', 'South America')
env = gentry_env[gentry_env$continent %in% gd_contintents, ]
# filter site table down to only samples in cts
env = env[env$abbreviation %in% unique(cts$site_code), ]
# filter cts table down to same set of sites
row_indices = cts$site_code %in% env$abbreviation
cts = cts[row_indices, ]

gentry_clean = merge(cts, env, by.x='site_code', by.y='abbreviation',
                     all=T)

map('world')
points(gentry_clean$lon, gentry_clean$lat, pch=1, col='red')

write.csv(gentry_clean, file='./data/filtered_data/gentry.csv',
          row.names=F)

# analysis ----------------------------------------------------------------
dat = list()
dat$scbi = read.csv('./data/filtered_data/SCBI_initial_woody_stem_census_2012_filtered.csv')
head(dat$scbi)

#dat$huss = read.csv('./data/filtered_data/00_HUSS_TREE_2013_filtered.csv')
#head(dat$huss)

dat$bci = read.csv('./data/filtered_data/bci_census7_filtered.csv')
head(dat$bci)

#########################################

plot_size = c(40, 40)

domain = c(0, 400, 0, 640)
comms = list()
comms$scbi = stems2comm(dat[[1]]$Latin, dat[[1]][ , c('gx', 'gy')], plot_size, domain)

domain = c(180, 820, 50, 450)
comms$bci = stems2comm(dat[[2]]$Latin, dat[[2]][ , c('gx', 'gy')], 
                       plot_size, domain) 

# append data together into one matrix for analysis
Spool_diff = ncol(comms$bci$comm) - ncol(comms$scbi$comm)
comms$scbi$comm = cbind(comms$scbi$comm, 
                        matrix(0, ncol=Spool_diff, 
                               nrow=nrow(comms$scbi$comm)))

comm = rbind(comms$scbi$comm, comms$bci$comm)
spat = rbind(comms$scbi$spat, comms$bci$spat)
env = data.frame(group=as.factor(rep(c('scbi', 'bci'), each=160)))

comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))


# not sure what is broken but
names(comm$env) = 'groups'


tst = get_delta_stats(comm, 'groups', ref_group='bci', type='discrete', 
                      log_scale=TRUE, inds = 10, nperm=10)

save(tst, file='./results/tst.Rdata')

######
## Jon's data
#########
# 1) invasion data
require(R.matlab)
dat_dir = paste('./data/', 'joninvade.mat', sep = '')

dat_matlab = readMat(dat_dir)
dat_plot = as.data.frame(matrix(NA, length(dat_matlab$x), 5))
dat_sp = as.data.frame(matrix(NA, length(dat_matlab$x), nrow(dat_matlab$comary) + 1))

dat_plot[, 1] = 1:nrow(dat_plot)
dat_plot[, 2] = as.vector(ifelse(unlist(dat_matlab$is.invaded) > 0, 'invaded', 'uninvaded'))
dat_plot[, 3] = unlist(dat_matlab$x)
dat_plot[, 4] = unlist(dat_matlab$y)
dat_plot[, 5] = rep(1, nrow(dat_plot))

dat_sp[, 1] = 1:nrow(dat_plot)
dat_sp[, 2:ncol(dat_sp)] = t(dat_matlab$comary)

# give reasonable names to env. data set:

comm = dat_sp
spat = dat_plot[,3:4]
env = dat_plot[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))


names(comm$env) = 'groups'


tst.inv = get_delta_stats(comm, 'groups', ref_group='uninvaded', type='discrete', log_scale=TRUE, inds=NULL, nperm=100)

save(tst.inv, file='./results/tst.inv.Rdata')

pdf("tst.inv.pdf", height = 5, width = 10)
plot.mobr(tst.inv)
dev.off()

pdf("tst.inv.box.pdf", height = 6, width = 8)
boxplot.comm(comm, "groups")
dev.off()

pdf("tst.inv.rarefy.pdf", height = 5, width = 10)
plot_rarefy(tst.inv)
dev.off()

#############
# 2) Morlaix #
#############

dat_dir = paste('./data/', 'morlaix.mat', sep = '')

dat_matlab = readMat(dat_dir)

dat_plot = as.data.frame(matrix(NA, 6, 5))
dat_sp = as.data.frame(matrix(NA, 6, nrow(dat_matlab$ab) + 1))

dat_plot[, 1] = 1:nrow(dat_plot)
dat_plot[, 2] = c(rep('before', 3), rep('after', 3))
dat_plot[, 3] = rep(1, nrow(dat_plot))
dat_plot[, 4] = rep(1, nrow(dat_plot))
dat_plot[, 5] = rep(1, nrow(dat_plot))

dat_sp[, 1] = 1:nrow(dat_plot)
dat_sp[, 2:ncol(dat_sp)] = t(dat_matlab$ab[, c(1:3, 6:8)])

comm = dat_sp
env = dat_plot[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env))
names(comm$env) = 'groups'

tst.mor = get_delta_stats(comm, 'groups', ref_group='before', type='discrete', log_scale=TRUE, inds=NULL, nperm=100)
save(tst.mor, file='./results/tst.mor.Rdata')

pdf("tst.mor.pdf", height = 5, width = 10)
plot.mobr(tst.mor)
dev.off()

pdf("tst.mor.box.pdf", height = 6, width = 8)
boxplot.comm(comm, "groups")
dev.off()

pdf("tst.mor.rarefy.pdf", height = 5, width = 10)
plot_rarefy(tst.mor)
dev.off()

#####
# 3) Jon's fire data:

dat_unburned = read.csv("./data/fire_data_unburned.csv")
head(dat_unburned)
colnames(dat_unburned)
dim(dat_unburned)

dat_burned = read.csv("./data/fire_data_burned.csv")
head(dat_burned)
colnames(dat_burned)
dim(dat_burned)

names(dat_burned)[2] <- "Treatment"
dat_unburned$Plot <- toupper(dat_unburned$Plot)
dat_burned$Plot <- tolower(dat_burned$Plot)

# rbind rows

library(dplyr)

dat.all <- bind_rows(dat_unburned, dat_burned)
str(dat.all)
names(dat.all)

dat.all.2 <- dat.all[,-c(23)]
tail(dat.all.2)

dat.all.2[is.na(dat.all.2)] <- 0
dat.all.2$Treatment <- as.character(dat.all.2$Treatment)
dat.all.2$Treatment[25:48] <- "burned"

## get coordinates:
dat_burned_xy= read.csv("./data/Fire_lat_longs.csv")
names(dat_burned_xy)[1] <- "Plot"

# change first letter into capital
dat_burned_xy$Plot <- as.character(dat_burned_xy$Plot)

dat_burned_xy$Plot[1:24] <- tolower(dat_burned_xy$Plot)[1:24]
dat_burned_xy$Plot[25:48] <- toupper(dat_burned_xy$Plot)[25:48]

match(sort(dat.all.2$Plot), sort(dat_burned_xy$Plot))

sort(dat.all.2$Plot)
sort(dat_burned_xy$Plot)

## which coordinates are potentially strange

dat <- inner_join(dat.all.2[,1:2], dat_burned_xy)

comm = dat.all.2[,3:23]
spat = dat_burned_xy[,2:3]
env = dat.all.2[,2]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))

names(comm$env) = 'groups'

tst.fire = get_delta_stats(comm, 'groups', ref_group='unburned', type='discrete', log_scale=TRUE, inds=NULL, nperm=100)
save(tst.fire, file='./results/tst.fire.Rdata')

pdf("tst.fire.pdf", height = 5, width = 10)
plot.mobr(tst.fire)
dev.off()

pdf("tst.fire.box.pdf", height = 6, width = 8)
boxplot.comm(comm, "groups")
dev.off()

pdf("tst.fire.rarefy.pdf", height = 5, width = 10)
plot_rarefy(tst.fire)
dev.off()


#####
# 4) jon coffee
######

dat_coffee = read.csv("./data/coffee_comm.csv")
head(dat_coffee)
colnames(dat_coffee)
dim(dat_coffee)

dat_xy = read.csv("./data/coffee_xy.csv")
dat_xy$Treatment <- rep(c("Natural", "Shaded"), each = 3)
head(dat_xy)
colnames(dat_xy)
dim(dat_xy)


comm = dat_coffee[, 2:22]
spat = dat_xy[,2:3]
env = dat_xy[,4]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))

names(comm$env) = 'groups'

tst.coffee = get_delta_stats(comm, 'groups', ref_group='Natural', type='discrete', log_scale=TRUE, inds=NULL, nperm=1000)

save(tst.coffee, file='./results/tst.coffee.Rdata')

pdf("tst.coffee.pdf", height = 5, width = 10)
plot.mobr(tst.coffee)
dev.off()

pdf("tst.coffee.box.pdf", height = 6, width = 8)
boxplot.comm(comm, "groups")
dev.off()

pdf("tst.coffee.rarefy.pdf", height = 5, width = 10)
plot_rarefy(tst.coffee)
dev.off()

###
# 5) Cattle_tank
####

dat_cattle_high= read.csv("./data/Cattle_tank_high.csv")
head(dat_cattle_high)
colnames(dat_cattle_high)
dim(dat_cattle_high)
dat_cattle_high[is.na(dat_cattle_high)] <- 0

rownames(dat_cattle_high) <- dat_cattle_high[,1]
dat_cattle_high <- dat_cattle_high[,-1]
dat_cattle_high <- t(dat_cattle_high)


dat_cattle_low= read.csv("./data/Cattle_tank_low.csv")
head(dat_cattle_low)
colnames(dat_cattle_low)
dim(dat_cattle_low)
dat_cattle_low[is.na(dat_cattle_low)] <- 0

rownames(dat_cattle_low) <- dat_cattle_low[,1]
dat_cattle_low <- dat_cattle_low[,-1]
dat_cattle_low <- t(dat_cattle_low)

identical(colnames(dat_cattle_low), colnames(dat_cattle_high))

dat_cattle <- rbind(dat_cattle_low, dat_cattle_high)
rownames(dat_cattle) <- paste0(dat_cattle_xy[,2], dat_cattle_xy[,1])

## coordinates:
dat_cattle_xy= read.csv("./data/Cattle_tank_xy.csv")
rownames(dat_cattle_xy) <- rownames(dat_cattle)

comm = dat_cattle
spat = dat_cattle_xy[,3:4]
env = dat_cattle_xy[,1]

## make community matrix:
comm = make_comm_obj(comm, data.frame(group=env, x=spat[,1], y=spat[,2]))
names(comm$env) = 'groups'

tst.cattle = get_delta_stats(comm, 'groups', ref_group='Low', type='discrete', log_scale=TRUE, inds=NULL, nperm=100)
save(tst.cattle, file='./results/tst.cattle.Rdata')

pdf("tst.cattle.pdf", height = 5, width = 10)
plot.mobr(tst.cattle)
dev.off()

pdf("tst.cattle.box.pdf", height = 6, width = 8)
boxplot.comm(comm, "groups")
dev.off()

pdf("tst.cattle.rarefy.pdf", height = 5, width = 10)
plot_rarefy(tst.cattle)
dev.off()
