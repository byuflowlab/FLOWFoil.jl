
#"""
#"""
#function initialize_wake()

#    # initialize the wake panel start point and direction
#    # wakestart = nodes[1]
#    # wakedir = [1.0; 0.0]

#    #TODO LATER (probably in different function)
#    # update the wake panel starting location to be the midpoint of the gap panel.
#    #        wakestart = (nodes[end] .+ nodes[end - 1]) / 2.0

#    # and update the wake panel initial direction to be the normal of that panel
#    #       wakedir = FLOWFoil.get_normal(nodes[end - 1], nodes[end])
#    #TODO: probably put this elsewhere
#    # update wake panel direction to be bisection of trailing edge panel vectors
#    # get vector along first panel
#    #        a1 = nodes[1] - nodes[2]

#    # get vector along second panel
#    #       an = nodes[end] - nodes[end - 1]

#    # calculate vector that bisects the first and last panel vectors
#    #      bisector = a1 * LinearAlgebra.norm(an) + an * LinearAlgebra.norm(a1)

#    # normalize to get the unit vector
#    #     wakedir = bisector / LinearAlgebra.norm(bisector)

#    #TODO: There is something in the method about a trailing half panel, find out what that means and if you should remove the final wake node or not.
#    # get initial wake geometry: equidistant panels starting at wakestart point and extending the input percentage of the airfoil chord in the calculated direction.
#    # wake_nodes = [
#    #    wakestart .+ x .* wakedir for
#    #     x in range(0.0; stop=wakelength * chordlength, length=numwake)
#    # ]

#    # get wake midpoints as well
#    #wake_midpoints = [
#    #    [(wake_nodes[i + 1][1] + wake_nodes[i][1]) / 2.0 (
#    #        wake_nodes[i + 1][2] + wake_nodes[i][2]
#    #    ) / 2.0] for i in 1:(numwake - 1)
#    #]

#end
