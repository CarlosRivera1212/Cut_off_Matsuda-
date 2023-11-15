# Define a function f_vn that determines an observation window in two coordinates based on groups.
# Arguments:
#   x: Numeric vector representing the x-coordinate.
#   y: Numeric vector representing the y-coordinate.
#   g: Factor vector specifying groups ('IR' and 'NIR').

f_vn = function(x, y, g) {
  
  # Find the maximum x-coordinate for elements where the group is 'IR'.
  px_mn = max(x[g == 'IR'])
  
  # Find the minimum x-coordinate for elements where the group is 'NIR'.
  px_mx = min(x[g == 'NIR'])
  
  # Calculate the midpoint of the x-coordinate range.
  px_m = (px_mn + px_mx) / 2
  
  # Find the minimum y-coordinate for each group ('IR' and 'NIR').
  pre_y_mn = rev(tapply(y, g, min))
  
  # Find the maximum y-coordinate for each group ('IR' and 'NIR').
  pre_y_mx = tapply(y, g, max)
  
  # Identify the group with the smallest y-coordinate range and its corresponding minimum and maximum y-coordinates.
  py_mn = pre_y_mn[which.min(pre_y_mx - pre_y_mn)]
  py_mx = pre_y_mx[which.min(pre_y_mx - pre_y_mn)]
  
  # Calculate the midpoint of the y-coordinate range.
  py_m = (py_mn + py_mx) / 2
  
  # Return a list containing the median x and y coordinates and the interval values.
  return(list(
    med = c(px_m, py_m),
    int = c(px_mn, px_mx, py_mn, py_mx)
  ))
}
