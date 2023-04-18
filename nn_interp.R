# Toy example of `Besag` spatial interpolation for spatial point data -----
# Tommy Irons, 04/2023, PhD @ University of Exeter.

# This script simulates normally distributed data as if they were point data
# observed at select locations in a spatial field. The object of the script
# is to create a function that can run the data through the analysis whilst
# allowing for parameter changes that will affect the final result.


# The point-data pattern will be triangulated, using a planar straight line
# graph. This is all mostly done using the packages RTriangle and spatstat- they're great.
# The model used is the (improper) Besag spatial model. Here, each point is
# given an expectation of the average of it's neighbours, alongside some
# Gaussian variance. This latent Gaussian model can be shown to be related
# to a Gauss-Markov Random Field, for which efficient computational inference
# methods already exist, as in this here script. The package Matrix is great
# for all of this computational linear algebra, and sparse matrices.

library(RTriangle) # Delaunay triangulation of point-data patterns
library(Matrix)    # Sparse matrices
library(tidyverse) # Wrangling
library(ggvoronoi) # Voronoi tesselation in ggplots

# Simulating the data to motivate the function ----------------------------
set.seed(123)

# Observation locations
coords1 = crossing(x=seq(0,5,len=10),y=seq(0,5,len=10)) %>%
  mutate(x = x + rnorm(n()),
         y = y + rnorm(n()))
# plot(coords1$x, coords1$y)

# Some basic trend (this could do with improvements, e.g Re(fft()))
orig = seq(0,11,len=100)

# Dataframe of coordinates and "truth" data
data = bind_cols(coords1, orig)
names(data) = c("x","y","orig")
#summary(data)

# Plotting original data ("truth")
#ggplot(data) + geom_point(aes(x=x,y=y,colour=orig))
#ggplot(data) + geom_voronoi(aes(x=x,y=y,fill=orig))

# Adding noise to "truth" to simulate noisy observations
data = data %>%
  mutate(z = orig + rnorm(n()))

#Plotting this
orig_mean_plot = ggplot(data) +
  geom_voronoi(aes(x=x,y=y,fill=z)) +
  geom_point(aes(x=x,y=y), col="red", data = data) +
  labs(title = "Observations") + theme_bw()
orig_mean_plot
# A little noisier. Motivation for smoothing :)


# Define coordinates as a Matrix for construction of planar straight line graph
coords = as.matrix(data[,c("x","y")])

# Standard deviations (just same trend scaled by 1/4)
var_y = seq(0,11/4, len=101)[2:101] # non-zero !

# Variance-covariance matrix
Sigma_y = sparseMatrix(
  i = 1:nrow(data),
  j = 1:nrow(data),
  x = var_y)

# That's all we need.


nn_interp = function(data,                      # 1-column matrix
                     int_col,                   # Number > 0
                      coords,                   # 2-column matrix
                      sigma_mat,                # NxN matrix for N observations
                      steiners = nrow(coords),  # Integer
                      tri_param = c(0.1, 25),   # 2-vector (area, angle)
                      tau = 1){                 # Number > 0
  
  # Things you must provide:
  #   - The vector of observations                                  (data)
  #   - The coordinate pairs of the observations                    (coords)
  #   - The covariance matrix of the observations                   (sigma_mat)
  
  # Things you have control over:
  #   - The number of additional points added                       (steiners)
  #   - The triangulation parameters (angle, area etc.)             (tri_param)
  #   - Tau, the smoothing parameter for the latent field           (tau)
  
  
  
  
  
  ## PS: THIS SCRIPT INCLUDES ANALYSIS OF IRIGNAL DATA WHICH IS NOT INCLUDED IN INPUT ##
  cat("Boring stuff...", "\n")
  
  # Libraries #
  library(tidyverse); library(ggvoronoi); library(gridExtra) # Data wrangling
  library(Matrix) # Efficient matrix computation
  library(RTriangle) # Very good triangulation of point-data package
  
  
  # Data # 
  # Plotting data as voronoi tesselation of the points
  orig_mean_plot = ggplot(data) +
    geom_voronoi(aes(x=x,y=y,fill=orig)) +
    geom_point(aes(x=x,y=y), col="red", data = data) +
    labs(title = "Original Simulations") + theme_bw()
  
  data_mean_plot = ggplot(data) +
    geom_voronoi(aes(x=x,y=y,fill=z)) +
    geom_point(aes(x=x,y=y), col="red", data = data) +
    labs(title = "`Noisy` Observations") + theme_bw()
  
  cat("Initiating triangulation.", "\n")
  # Creating planar straight line graph (object to store point data)
  point_graph = pslg(coords)
  
  # Getting triangulation parameters from input
  tri_area = tri_param[1] # Maximum triangle area
  tri_angle = tri_param[2] # Minimum triangle area
  tri_steiners = steiners # Maximum number of added (Steiner) points (latent field)
  
  # Creating triangulation object for the specified parameters
  tri_obj = RTriangle::triangulate(
    p = point_graph,
    a = tri_area,
    q = tri_angle,
    S = tri_steiners  )
  
  cat("Triangulation complete.", "\n")
  # Number of points in final triangulation (total points in latent field)
  N_p = nrow(tri_obj$P)
  
  cat("Constructing matrices...", "\n")
  # Create an object that stores the coordinates of the points in tri_obj
  A_tib = tibble(
    ind = 1:N_p,
    x = tri_obj$P[,1],
    y = tri_obj$P[,2])
  
  # Append A_tib to data, assigning NA at the new points
  df = full_join(A_tib, data)
  
  # Getting indices of the observed points
  inds = which(!is.na(df[,int_col]))
  
  # Creating A matrix to pull observations from latent field
  A = sparseMatrix(i = inds, # Rows of non-zero entries
                   j = inds, # Columns of non-zero entries
                   x = rep(1, length(inds)) )
  
  # Arranging neighbours in useful way for construction of adjacency matrix
  nb_list = data.frame(p1 = tri_obj$E[,1], p2 = tri_obj$E[,2]) %>%
    transform(p1 = pmin(p1, p2),
              p2 = pmax(p1, p2))
  
  # Constructing the adjacency matrix
  adj_mat = sparseMatrix(
    i = nb_list$p1, # `Reference point`
    j = nb_list$p2, # `Neighbour`
    x = 1,
    symmetric = TRUE)
  
  if (sum(dim(adj_mat))/2 == N_p){
    cat("Adjacency matrix looks okay at first glance...", "\n")
  } else {
    stop("Uh oh... adjacency matrix does not look good", "\n")
  }
  
  # Getting number of nearest neighbours
  adj_mat_sums = rowSums(adj_mat)
  
  # Define dataframe that tells us some stuff about our neighbour pairings
  nb_count = data.frame(table(tri_obj$E)) %>% # Frequency table of pairings
    rename(idx = Var1) %>% # index
    mutate(idx = as.numeric(idx) )
  
  # Row-wise coefficients for a matrix
  nb_vec = c(-1/nb_count$Freq)
  
  # Cholesky decomposition of posterior precision is this matrix
  C = nb_vec * adj_mat + diag(x=1, nrow = N_p, ncol = N_p)
  
  
  cat("Generating precision matrix for data...", "\n")
  Q_y = solve(sigma_mat)
  
  cat("Generating precision matrix for latent spatial prior...", "\n")
  # Var[X] = (tau * C'C)^{-1}
  tau = tau                                 # Hyperparameter
  Q_x = tau * crossprod(C,C)
  Var_X = solve(Q_x)
  
  z = as.matrix(data[,int_col])
  
  
  cat("Canonical mean...", "\n")
  b = crossprod(A, Q_y%*%z)
  
  
  cat("Canonical variance...", "\n")
  Q_post = Q_x + crossprod(A, Q_y%*%A)
  
  
  cat("Conditional expectation...", "\n")
  mean = solve(Q_post) %*% b
  
  
  cat("Conditional variance...", "\n")
  variance = solve(Q_post)
  
  # Storing conditional expectation and variance of LGM (given the data)
  post_mean = mean[,1]
  post_var = diag(variance)
  
  
  cat("Generating plots...", "\n")
  post_dat = df %>%
    mutate(post_mean = post_mean,
           post_var = post_var)
  
  post_mean_plot = ggplot(post_dat) +
    geom_voronoi(aes(x=x,y=y,fill=post_mean)) +
    labs(title = "Interpolation - mean") +
    theme_bw()
  
  post_var_plot = ggplot(post_dat) +
    geom_voronoi(aes(x=x,y=y,fill=post_var)) +
    labs(title = "Interpolation - variance") +
    theme_bw()
  
  
  cat("Here are some pretty plots!", "\n")
  grid.arrange(orig_mean_plot, data_mean_plot,
               post_mean_plot, post_var_plot,
               ncol = 2)
  
  return(post_mean)
}


# Basic function test first
test = nn_interp(data = data, coords=coords, int_col = 4,
          sigma_mat = Sigma_y,
          steiners = 0,
          tri_param = c(0.15, 25),
          tau = 10)
# Takes less than a second





