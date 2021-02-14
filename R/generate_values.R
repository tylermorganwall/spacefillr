#'@title Generate Halton Value (Randomly Initialized)
#'
#'@description Generate a single value from a seeded Halton set.
#'
#'Note: This is much slower than generating the entire set ahead of time.
#'
#'@param i The element of the sequence to extract.
#'@param dim The dimension of the sequence to extract.
#'@param seed Default `0`. The random seed.
#'@return A single numeric value representing the `i`th element in the `dim` dimension.
#'
#'@export
#'@examples
#'#Generate a 3D sample:
#'point3d = c(generate_halton_random_single(10, dim = 1),
#'            generate_halton_random_single(10, dim = 2),
#'            generate_halton_random_single(10, dim = 3))
#'point3d
#'
#'#Change the random seed:
#'#'#Generate a 3D sample
#'point3d_2 = c(generate_halton_random_single(10, dim = 1, seed = 10),
#'              generate_halton_random_single(10, dim = 2, seed = 10),
#'              generate_halton_random_single(10, dim = 3, seed = 10))
#'point3d_2
generate_halton_random_single = function(i, dim, seed = 0) {
  return(rcpp_generate_halton_random_single(i,dim,seed))
}

#'@title Generate Halton Value (Faure Initialized)
#'
#'@description Generate a single value from a seeded Halton set, initialized with a Faure sequence.
#'
#'Note: This is much slower than generating the entire set ahead of time.
#'
#'@param i The element of the sequence to extract.
#'@param dim The dimension of the sequence to extract.
#'@return A single numeric value representing the `i`th element in the `dim` dimension.
#'
#'@export
#'@examples
#'#Generate a 3D sample:
#'point3d = c(generate_halton_faure_single(10, dim = 1),
#'            generate_halton_faure_single(10, dim = 2),
#'            generate_halton_faure_single(10, dim = 3))
#'point3d
generate_halton_faure_single = function(i, dim) {
  return(rcpp_generate_halton_faure_single(i,dim))
}

#'@title Generate Halton Set (Randomly Initialized)
#'
#'@description Generate a set of values from a seeded Halton set.
#'
#'@param n The number of values (per dimension) to extract.
#'@param dim The number of dimensions of the sequence.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `dim` matrix listing all the
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_halton_random_set(n=1000, dim=2)
#'plot(points2d)
#'
#'#Change the seed and extract a separate pair of dimensions
#'points2d = generate_halton_random_set(n=1000, dim=10,seed=2)
#'plot(points2d[,5:6])
#'
#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = matrix(generate_halton_random_set(10000,dim=2),ncol=2)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_halton_random_set = function(n, dim, seed = 0) {
  vals = unlist(rcpp_generate_halton_random_set(n,dim,seed))
  return(matrix(vals, nrow=n,ncol=dim))
}

#'@title Generate Halton Set (Faure Initialized)
#'
#'@description Generate a set of values from a Faure Halton set.
#'
#'@param n The number of values (per dimension) to extract.
#'@param dim The number of dimensions of the sequence.
#'@return An `n` x `dim` matrix listing all the
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_halton_random_set(n=1000, dim=2)
#'plot(points2d)
#'
#'#Extract a separate pair of dimensions
#'points2d = generate_halton_random_set(n=1000, dim=10)
#'plot(points2d[,5:6])
#'
#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = matrix(generate_halton_faure_set(10000,dim=2),ncol=2)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_halton_faure_set = function(n, dim) {
  vals = unlist(rcpp_generate_halton_faure_set(n,dim))
  return(matrix(vals, nrow=n,ncol=dim))
}

#'@title Generate Sobol Set
#'
#'@description Generate a set of values from a Sobol set.
#'
#'@param n The number of values (per dimension) to extract.
#'@param dim The number of dimensions of the sequence.
#'@param seed Default `0`. The random seed.
#'@return A single numeric value representing the `i`th element in the `dim` dimension.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_sobol_set(n=1000, dim = 2)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_sobol_set(n=1500, dim = 2)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = matrix(generate_sobol_set(10000,dim=2),ncol=2)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_sobol_set = function(n, dim, seed = 0) {
  vals = unlist(rcpp_generate_sobol_set(n, dim, seed))
  return(matrix(vals,ncol=2))
}

#'@title Generate Owen-scrambled Sobol Set
#'
#'@description Generate a set of values from an Owen-scrambled Sobol set.
#'
#'@param n The number of values (per dimension) to extract.
#'@param dim The number of dimensions of the sequence.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `dim` matrix with all the calculated values from the set.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_sobol_owen_set(n=1000, dim = 2)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_sobol_owen_set(n=1500, dim = 2)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = matrix(generate_sobol_owen_set(10000,dim=2),ncol=2)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_sobol_owen_set = function(n, dim, seed = 0) {
  vals = unlist(rcpp_generate_sobol_owen_set(n, dim, seed))
  return(matrix(vals,ncol=2))
}

#'@title Generate Owen-scrambled Sobol Value (fast approximate method)
#'
#'@description Generate a set of values from an Owen-scrambled Sobol set using an
#'approximate hashing method.
#'
#'@param n The number of values (per dimension) to extract.
#'@param dim The number of dimensions of the sequence.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `dim` matrix with all the calculated values from the set.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_sobol_owen_fast_set(n=1000, dim = 2)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_sobol_owen_fast_set(n=1500, dim = 2)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = matrix(generate_sobol_owen_fast_set(10000,dim=2),ncol=2)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_sobol_owen_fast_set = function(n, dim, seed = 0) {
  vals = unlist(rcpp_generate_sobol_owen_fast_set(n, dim, seed))
  return(matrix(vals,ncol=2))
}


#'@title Generate 2D Progressive Jittered Set
#'
#'@description Generate a set of values from a Progressive Jittered set.
#'
#'@param n The number of 2D values to extract.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `2` matrix with all the calculated values from the set.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_pj_set(n=1000)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_pj_set(n=1500)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a new set by changing the seed
#'points2d = generate_pj_set(n=1500,seed=10)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = generate_pj_set(10000)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_pj_set = function(n, seed = 0) {
  vals = unlist(rcpp_generate_pj_set(n, seed))
  return(matrix(vals, nrow=n,ncol=2,byrow = TRUE))
}

#'@title Generate 2D Progressive Multi-Jittered Set
#'
#'@description Generate a set of values from a Progressive Multi-Jittered set.
#'
#'@param n The number of 2D values to extract.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `2` matrix with all the calculated values from the set.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_pmj_set(n=1000)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_pmj_set(n=1500)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a new set by changing the seed
#'points2d = generate_pmj_set(n=1500,seed=10)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = generate_pj_set(10000)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_pmj_set = function(n, seed = 0) {
  vals = unlist(rcpp_generate_pmj_set(n, seed))
  return(matrix(vals, nrow=n,ncol=2,byrow = TRUE))
}

#'@title Generate 2D Progressive Multi-Jittered (with blue noise) Set
#'
#'@description Generate a set of values from a Progressive Multi-Jittered (with blue noise) set.
#'
#'@param n The number of 2D values to extract.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `2` matrix with all the calculated values from the set.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_pmjbn_set(n=1000)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_pmjbn_set(n=1500)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a new set by changing the seed
#'points2d = generate_pmjbn_set(n=1500,seed=10)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = generate_pmjbn_set(10000)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_pmjbn_set = function(n, seed = 0) {
  vals = unlist(rcpp_generate_pmjbn_set(n, seed))
  return(matrix(vals, nrow=n,ncol=2,byrow=TRUE))
}

#'@title Generate 2D Progressive Multi-Jittered (0, 2) Set
#'
#'@description Generate a set of values from a Progressive Multi-Jittered (0, 2) set.
#'
#'@param n The number of 2D values to extract.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `2` matrix with all the calculated values from the set.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_pmj02_set(n=1000)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_pmj02_set(n=1500)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a new set by changing the seed
#'points2d = generate_pmj02_set(n=1500,seed=10)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = generate_pmj02_set(10000)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_pmj02_set = function(n, seed = 0) {
  vals = unlist(rcpp_generate_pmj02_set(n, seed))
  return(matrix(vals, nrow=n,ncol=2,byrow=TRUE))
}

#'@title Generate 2D Progressive Multi-Jittered (0, 2) (with blue noise) Set
#'
#'@description Generate a set of values from a Progressive Multi-Jittered (0, 2) (with blue noise) set.
#'
#'@param n The number of 2D values to extract.
#'@param seed Default `0`. The random seed.
#'@return An `n` x `2` matrix with all the calculated values from the set.
#'
#'@export
#'@examples
#'#Generate a 2D sample:
#'points2d = generate_pmj02bn_set(n=1000)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a longer sequence of values from that set
#'points2d = generate_pmj02bn_set(n=1500)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Generate a new set by changing the seed
#'points2d = generate_pmj02bn_set(n=1500,seed=10)
#'plot(points2d, xlim=c(0,1),ylim=c(0,1))
#'
#'#Integrate the value of pi by counting the number of randomly generated points that fall
#'#within the unit circle.
#'pointset = generate_pmj02bn_set(10000)
#'
#'pi_estimate = 4*sum(pointset[,1] * pointset[,1] + pointset[,2] * pointset[,2] < 1)/10000
#'pi_estimate
generate_pmj02bn_set = function(n, seed = 0) {
  vals = unlist(rcpp_generate_pmj02bn_set(n,seed))
  return(matrix(vals, nrow=n,ncol=2,byrow=TRUE))
}

