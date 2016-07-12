source('./mobr/R/mobr.R')

gentry = read.csv('./data/filtered_data/gentry.csv')

site_line = paste(gentry$site_id, gentry$line, sep='_')
gentry_comm = tapply(gentry$count, 
                     list(site_line, gentry$species_id),
                     sum)
gentry_comm = ifelse(is.na(gentry_comm), 0, gentry_comm)

row_indices = match(row.names(gentry_comm), site_line)
gentry_env = gentry[row_indices, 
                    c('site_id', 'line', 'country', 'continent',
                      'lon', 'lat', 'min_elev', 'max_elev',
                      'precip')]
row.names(gentry_env) = paste(gentry_env$site_id, gentry_env$line, sep='_')

gentry_comm = make_comm_obj(gentry_comm, gentry_env) 
                            
gentry_tst = get_delta_stats(gentry_comm, 'lat', type='continuous',
                             log_scale=T, inds=10, nperm=10) 
                              