# Define a function f_po that finds the values of "p" that minimize the distance
# Arguments:
#   x: Numeric vector representing the x-coordinate.
#   y: Numeric vector representing the y-coordinate.
#   g: Grouping variable.
#   mrx: Mean value of x.
#   mry: Mean value of y.
#   mvx: Target mean value for x.
#   mvy: Target mean value for y.

f_po = function(x, y, g, mrx, mry, mvx, mvy) {
  
  # Calculate the squared Spearman correlation coefficient between x and y.
  r = cor(x, y, method = 'spearman') ** 2
  
  # Calculate the pairwise distances between the means of x and y using f_mm function.
  d = as.matrix(dist(cbind(f_mm(x), f_mm(y))))
  di = 1 / d
  di[is.infinite(di)] = 0
  
  # Initialize an empty vector to store differences.
  dif = NULL
  
  # Loop over potential values of "p".
  for (i in pot) {
    dip = di ** i
    W = dip / apply(dip, 1, sum)
    W[is.na(W)] = 0
    
    # Weighted mean of x considering correlation.
    xm = mean((W %*% x) * r + (1 - r) * x)
    
    # Weighted mean of y considering correlation.
    ym = mean((W %*% y) * r + (1 - r) * y)
    
    # Calculate the Euclidean distance between (xm, ym) and (mvx, mvy).
    dif_i = sqrt((xm - mvx) ** 2 + (ym - mvy) ** 2)
    
    # Append the difference to the vector.
    dif = append(dif, dif_i)
  }
  
  # Select the "p" value that minimizes the differences.
  pot_sel = pot[which.min(dif)]
  dip = di ** pot_sel
  W = dip / apply(dip, 1, sum)
  W[is.na(W)] = 0
  
  # Calculate the weighted means of x and y using the selected "p".
  Wx = (W %*% x) * r + (1 - r) * x
  Wy = (W %*% y) * r + (1 - r) * y
  
  # Return a vector containing the mean values.
  return(c(mwx = mean(Wx), mwy = mean(Wy)))
}
