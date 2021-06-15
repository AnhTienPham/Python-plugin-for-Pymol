from pymol import cmd, stored

# A function to identify surface area of a protein based on the difference between selections
# and objects

def identifySA(pdb, chain, start, end):
    # Create lists to store the selected, object, difference of area and list of residue position
    selArea = []
    objArea = []
    diffArea = []
    listPos = []
    start = int(start)
    end = int(end)

    # Fetch the protein
    cmd.fetch(pdb)

    # Create selection and object
    cmd.select("sel_nanobody", pdb + " and " + chain)
    cmd.create("obj_nanobody", pdb + " and " + chain)

    # Create selection for all residues
    for x in range(start, end):
        cmd.select("sele" + str(x), "sel_nanobody and resi " + str(x) + " and not name c+n+o")
    # Find the surface areas of all the selections
    for x in range(start, end):
        selArea.append(cmd.get_area("sele" + str(x)))
    # Create objects for all residues
    for x in range(start, end):
        cmd.select("obj" + str(x), "obj_nanobody and resi " + str(x) + " and not name c+n+o")
    # Find the surface areas of all the objects
    for x in range(start, end):
        objArea.append(cmd.get_area("obj" + str(x)))
    # Find the difference in surface areas of all the residues
    for x in range(0, end - 1):
        diffArea.append(selArea[x] - objArea[x])
    # If there is a difference in the area, then print the position of the residues
    for x in range(0, len(diffArea)):
        if diffArea[x] < -0.5:
            listPos.append(x + 1)
    print(listPos)

cmd.extend( "identifySA", identifySA );