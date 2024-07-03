"""
Code to generate a mesh of a domain with fractures
"""
# -----------------------
#   Importing modules
# -----------------------
from TPZMeshModeling import TPZMeshModeling
import gmsh

# -----------------------
#       Functions
# -----------------------
def DispX0(x0:tuple[float], disp:tuple[float])->list[float]:
    """
    Return the point x0 displaced by the vector disp
    @carlos : what does this statement mean?
    """
    return [X + d for X, d in zip(x0, disp)]

# -----------------------
#   Main function
# -----------------------
def main()->None:
    # ATTENTION: remeber to always create the mesh during the show model 
    # since now the model is shown three times, it is possible that the mesh is not created in one of them. 
    # Just type 2 on your keyboard to be sure that the mesh is created

    # initial coordinate
    x0 = [0., 0., 0.]

    # length values
    d1:float = 1 # domain 1 length 
    d2:float = 1 # domain 1 height
    d3:float = d2 / 2 # fracture height
    d4:float = d1 / 2 # fracture length

    # mesh size
    lc = 2e-1

    # show mesh bool 
    # change this to hide or show the mesh 
    show = True

    TPZMeshModeling.Begin()
    TPZMeshModeling.TurnOnLabels('point', 'line', 'surface') # turn on the labels of the points, lines and surfaces
    TPZMeshModeling.TurnOnRepresentation('surfaces') # turn on the geometrical representation of the surfaces

    # Creating the surface domain
    # point coordinates
    # @carlos : what is M1_coord and M2_coord?
    rectangle_points = [
        x0, # p1
        DispX0(x0, (d1, 0., 0.)), # p2
        DispX0(x0, (d1, d2, 0.)), # p3
        DispX0(x0, (0, d2, 0.)), # p4
        DispX0(x0, (0, d3, 0)) # p5
    ]
    # @carlos : please indicate on the drawing where these points are?
    p1, p2, p3, p4, p5= TPZMeshModeling.CreatePoints(rectangle_points, lc)

    # lines
    line_points = [
        [p1, p2], # l1 
        [p2, p3], # l2
        [p3, p4], # l3
        [p4, p5], # l4
        [p5, p1] # l5
        ]    

    # @carlos : please indicate where these lines are on the drawing?
    # I thought l1 , l2, etc where float values?
    l1, l2, l3, l4, l5 = TPZMeshModeling.CreateLines(line_points)

    # curve loops
    c1 = TPZMeshModeling.CreateCurveLoops([[l1, l2, l3, l4, l5]]) # domain 1 curve loop 

    # plane 
    plane1 = TPZMeshModeling.CreatePlanes([c1])[0] # domain 1 plane surface

    # fracture points
    frac_points = [
        (d4, d3, 0)
    ]

    fp2 = TPZMeshModeling.CreatePoints(frac_points, lc)[0]

    # fracture lines
    fl1 = TPZMeshModeling.CreateLines([[p5,fp2]])[0] # fracture 1 line

    TPZMeshModeling.Synchronize()

    gmsh.model.mesh.embed(1, [fl1], 2, plane1) # embed the fracture lines in the domain 1 plane

    TPZMeshModeling.Synchronize()

    if show:
        TPZMeshModeling.ShowModel() # showing the model without the physical groups

    # physical groups
    pg = [
        [(2, [plane1]), 1, "dom"],
        [(1, [l1]), 2, "BcB"],
        [(1, [l2]), 3, "BcR"],
        [(1, [l3]), 4, "BcT"],
        [(1, [l4,l5]), 5, "BcL"],
        [(1, [fl1]), 6, "Fracture"],
        [(0, [fp2]), 7, "FractureEdge"]
    ]

    TPZMeshModeling.CreatePhysicalGroup(pg) # creating the physical groups

    TPZMeshModeling.CreateMesh(2) # creating the mesh

    if show:
        # showing the model only with the physical group tags 
        TPZMeshModeling.TurnOffLabels('point')
        TPZMeshModeling.SetDescription("physical_tag")
        TPZMeshModeling.ShowModel()

        # showing the model only with the physical group names
        TPZMeshModeling.SetDescription("physical_name")
        TPZMeshModeling.ShowModel()

    TPZMeshModeling.WriteMeshFiles("SquareFrac", ".msh")

    TPZMeshModeling.End()


if __name__ == "__main__":
    main()