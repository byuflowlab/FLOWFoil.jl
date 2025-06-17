# """
#     surfacenormal(xloc, x, y)

# Finds x,y location and angle corresponding to surface of airfoil at specified
# x-location.
# """
# function surface_normal_angle(xloc, x, y)
#     idx = indmin(abs(x - xloc))
#     xloc = x[idx]
#     yloc = y[idx]
#     theta = atan2((y[idx + 1] - y[idx - 1]), ((x[idx + 1] - x[idx - 1])))
#     return xloc, yloc, theta
# end
