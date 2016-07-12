library(maps)
library(ecoretriever)
#library(mobr)
source('./mobr/R/mobr.R')
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

range(dat$gx)
range(dat$gy)


write.csv(dat, file=paste('./data/filtered_data/', fileprefix, '_filtered.csv', sep=''), 
          row.names=F)

# huss --------------------------------------------------------------------

fileprefix = '00_HUSS_TREE_2013'
dat = read.csv(paste('./data/', fileprefix, '.csv', sep=''), sep=';')

# fix typo
dat$SPECIES[dat$SPECIES == 'Bu'] = 'BU'

good_data = dat$PROBLEM != 1 &       # drop 536 indiv
  dat$SPECIES != '-1000'   # drop 2 indiv
dat = dat[good_data, ]

write.csv(dat, file=paste('./data/filtered_data/', fileprefix, '_filtered.csv', sep=''), 
          row.names=F)


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
line_cts = tapply(cts$line, list(as.character(cts$site_code)), 
                  function(x) length(unique(x)))
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
dat$bci = read.csv('./data/filtered_data/bci_census7_filtered.csv')
head(dat$bci)

plot_size = c(40, 40)
domain = c(0, 400, 0, 640)

comms = list()
comms$scbi = stems2comm(dat[[1]]$Latin, dat[[1]][ , c('gx', 'gy')], 
                        plot_size, domain) 
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
                      log_scale=TRUE, inds = 10, nperm=1000)

save(tst, file='./results/tst.Rdata')
