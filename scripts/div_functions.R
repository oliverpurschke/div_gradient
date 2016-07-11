

stems2comm = function(spids, stem_coords, plot_size, domain,
                      abu=NULL, rm_absent_sp=TRUE)
{ 
  ## Convert a mapped stem dataset into a community data matrix. 
  ## Output: 
  ## A community matrix where each row is a differnet pixel on a grid.  
  ## Arguments:
  ## spids: an integer or string specifying species identities
  ## stem_stem_coords : two column matrix (x,y) specifying the spatial coordinates of each stem
  ## plot_size : the x and the y size of each sub-plot.
  ## domain : specifies the spatial domain of the area:  (xmin, xmax, ymin, ymax)
  ## abu: abundance associated with each record, if NULL then it is set to 1
  ##      individual per record
  ## grainSuffix : if supplied the grain column will have this appended to it
  ##               so that it is clear what community this corresponds with
  ## rm_absent_sp: boolean that defaults to TRUE to remove species
  ##  who no longer occur in the site x species matrix after subsetting base
  ## on the defined spatial domain (i.e., argument 'domain' specifies a smaller
  ## area than spatial coordiantes are provided for)
  xdiff = abs(domain[2] - domain[1])
  ydiff = abs(domain[4] - domain[3])
  xlength = plot_size[1]
  ylength = plot_size[2]
  n_plots = (xdiff / xlength) * (ydiff / ylength)
  if (n_plots != round(n_plots)) 
    stop('number of plots that study site can be devided into is not an integer please chose a new plot_size')
  S = length(unique(spids))
  comm = matrix(NA, nrow=n_plots, ncol=S)
  spat = matrix(NA, nrow=n_plots, ncol=2)
  irow = 1
  xbreaks = seq(domain[1], domain[2], xlength)
  ybreaks = seq(domain[3], domain[4], ylength) 
  for (i in 1:(length(xbreaks) - 1)) {
    for (j in 1:(length(ybreaks) - 1)) {
      spat[irow, ] = c(xbreaks[i], ybreaks[j])
      if (i == length(xbreaks) - 1)
          xbreaks[i + 1] = xbreaks[i + 1] + 0.01
      if (j == length(ybreaks) - 1)
          ybreaks[j + 1] = ybreaks[j + 1] + 0.01
      in_plot =  xbreaks[i] <= stem_coords[ , 1] &
                 stem_coords[ , 1] < xbreaks[i + 1] & 
                 ybreaks[j] <= stem_coords[ , 2] &
                 stem_coords[ , 2] < ybreaks[j + 1]
      if (is.null(abu) ){
        comm[irow, ] = as.integer(table(c(spids[in_plot], 1:S)) - 1)
      }
      else {
        comm[irow, ] =  as.integer(table(c(unlist(mapply(
                                  rep, spids[in_plot], abu[in_plot])), 1:S)) - 1)
      }  
      irow = irow + 1 
    }
  }
  if (rm_absent_sp) {
    cols_to_rm = which(apply(comm, 2, function(x) all(x == '0')))
    if (length(cols_to_rm) > 0)
      comm = comm[ , -cols_to_rm]
  }
  params = data.frame(grain=prod(plot_size), extent=xdiff*ydiff, 
                      n_plots) 
  out = list(params=params, comm=comm, spat=spat)
  return(out)
}
