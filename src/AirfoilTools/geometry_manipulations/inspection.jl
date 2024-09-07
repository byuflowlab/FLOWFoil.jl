# """
#     surfacenormal(xloc, x, z)

# Finds x,z location and angle corresponding to surface of airfoil at specified
# x-location.
# """
# function surface_normal_angle(xloc, x, z)
#     idx = indmin(abs(x - xloc))
#     xloc = x[idx]
#     zloc = z[idx]
#     theta = atan2((z[idx + 1] - z[idx - 1]), ((x[idx + 1] - x[idx - 1])))
#     return xloc, zloc, theta
# end
