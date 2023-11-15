# Define a function f_mm that performs min-max normalization on a numeric vector x.
# Min-max normalization scales the values of x to a range between 0 and 1.
# Arguments:
#   x: Numeric vector to be normalized.

f_mm = function(x) {
  # Calculate the min-max normalization of x.
  normalized_x = (x - min(x)) / (max(x) - min(x))
  
  # Return the normalized vector.
  return(normalized_x)
}
