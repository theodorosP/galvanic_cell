from ase.io          import read, write
from dlePy.supercell import supercell, remove_atoms_outside
import numpy as np
import sys
from ase.io.pov import get_bondpairs, get_hydrogenbonds, set_high_bondorder_pairs
from ase.visualize import view
import matplotlib.cm as cm
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

def get_atom_setting( ):
    with open( '/shared/apps/VESTA-x86_64/elements.ini', 'r' ) as f:
        lines = f.readlines( )
    elements = { }
    for line in lines:
        tmp = line.split( )
        color = ( float( tmp[ 5 ] ), float( tmp[ 6 ] ), float( tmp[ 7 ] ) )
        rcov   = float( tmp[ 2 ] )
        rvdw   = float( tmp[ 3 ] )
        r   = float( tmp[ 4 ] )
        element = tmp[ 1 ]
        elements[ element ] = { }
        elements[ element ][ 'r'    ]  = r
        elements[ element ][ 'rvdw' ]  = rvdw
        elements[ element ][ 'rcov' ]  = rcov
        elements[ element ][ 'color' ] = color
    return elements


def get_figure_2( sys0, fout, rot = "-30x", w = 15, h = 15 ):
    system  = supercell( sys0, 1, 3, 1 )
    center = 0.5*sys0.cell[ 0 ] + 1.5*sys0.cell[ 1 ] + 0.5* sys0.cell[ 2 ]
    #center = 0.5*sys0.cell[ 0 ] + 0.5*sys0.cell[ 1 ] + 0.5* sys0.cell[ 2 ]
    #center[0] += sys0.cell[ 0, 0 ]
    center[2] -= 1.
    #center[2] += 3
    #center[1] += 10
    for i in range( 3 ):
        system.positions[ :, i ] -= center[ i ]
    
    #change color/radius
    elements = get_atom_setting( )
    elements[ 'S' ][ 'color' ] = ( 1, 0, 0 )
    elements[ 'O' ][ 'color' ] = ( 1, 0.7, 0.7 )
    elements[ 'H' ][ 'color' ] = ( 1, 1, 1 )
    elements[ 'C' ][ 'color' ] = ( 0.5, 0.5, 0.5 )
    elements[ 'He' ][ 'color' ] = ( 1, 0, 1 )
    elements[ 'S' ][ 'rcov' ] = elements[ 'O' ][ 'rcov' ]
    elements[ 'He' ][ 'rcov' ] = elements[ 'H' ][ 'rcov' ]
    elements[ 'F' ][ 'rcov' ] = elements[ 'O' ][ 'rcov' ]
    elements[ 'Ne' ][ 'rcov' ] = elements[ 'H' ][ 'rcov' ]
    elements[ 'F' ][ 'color' ] = elements[ 'O' ][ 'color' ]
    elements[ 'Ne' ][ 'color' ] = elements[ 'H' ][ 'color' ]
    elements[ 'Na' ][ 'rcov' ] = elements[ 'Na' ][ 'rcov' ] * 0.6

    del_list = []
    cutoff = 6.3
    for at in system:
        if at.symbol in [ "F", "Ne", 'C' ] and np.linalg.norm( at.position ) > cutoff:
            del_list.append( at.index )
    for at in sorted( del_list, reverse = True ):
        del system[ at ]

    radii = [ ]
    colors = [ ]
    for i, at in enumerate( system ):
        radius = elements[ at.symbol ][ 'rcov' ]
        if at.symbol in [ 'O', 'H', 'S', 'He', 'N' ]:
            radius = radius/3.
        radii.append( radius )
        color = elements[ at.symbol ][ 'color' ]
        alpha = 0
        if at.symbol in [ 'C', 'S', 'K' ]:
            alpha =1 
        factor = 0
        colors.append( ( factor + color[ 0 ], factor + color[ 1 ], factor + color[ 2 ], alpha ) )


    bond_pairs = get_bondpairs( system, radius = 0.70 )
    hydrogenbond1 = get_hydrogenbonds( system, atype1 = 'Na', atype2= [ 'O', 'S' ], radius = 4, rhbondrange = (2.3, 4.0) )

    hydrogenbond2 = get_hydrogenbonds( system, atype1 = 'H', atype2= [ 'O', 'S', 'N' ], radius = 5, rhbondrange = (1.3, 2.6) )
    hydrogenbond3 = get_hydrogenbonds( system, atype1 = 'He', atype2= [ 'O', 'S', 'N' ], radius = 5, rhbondrange = (1.3, 2.6) )
    #hydrogenbond2 = get_hydrogenbonds( system, atype1 = 'K', atype2= [ 'O', 'S' ], radius = 2, rhbondrange = (1.3, 3.0) )
    hydrogenbond  = {}
    for key in hydrogenbond1.keys():
        hydrogenbond[ key ] = hydrogenbond1[ key ]
    for key in hydrogenbond2.keys():
        hydrogenbond[ key ] = hydrogenbond2[ key ]
    for key in hydrogenbond3.keys():
        hydrogenbond[ key ] = hydrogenbond3[ key ]

    bond_pairs = set_high_bondorder_pairs(bond_pairs, high_bondorder_pairs=hydrogenbond)

    width  = w
    height = h
    bbox = (
           -width/2-1.5,
           -height/2-2.5,
           width/2-1.5,
           height/2-2.5,
           )
    #for rotx in [ 0, 15, 30, 45, 60, 75, 90 ]:
    #    rot = '90z,' + str( rotx ) + 'x' 
    write( fout+ '.pov', system, format = 'pov', run_povray = True,
           canvas_width = 1000,    # Set width, in pixel
           radii = radii,              # Set radius 
           bondatoms = bond_pairs, # Display bonds
           bbox  = bbox,
           colors = colors,        # Set colors
           celllinewidth = 0.0,
           rotation = rot,
           hydrogenbond = { 'ndots':5, 'color' : [0., 1, 0.], 'rdot':0.05 }
           )
    

if __name__ == "__main__":
    system = read( 'POSCAR4' )
    rot =  "270z,-80x"
    fout = "FS2"
    get_figure_2( system, fout, rot = "270z,-80x", w = 13, h = 12 )
