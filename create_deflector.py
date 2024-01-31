import gmsh
import numpy as np
import gmsh4mtc
import argparse

import os
import glob

"""
à écrire dans le terminal pour avoir le maillage :

python create_deflector.py --x 0.6 --H 0.5 --alpha 45
(les autres paramètres peuvent être appelés mais ce n'est pas nécessaire)

par défaut on crée un déflecteur arrondi, 
rajouter la commande : --simple pour avoir un déflecteur simple
"""

def get_args():
    
    parser = argparse.ArgumentParser(description='Description du script.')

    parser.add_argument('--x', type=float, required=True, help='position')
    parser.add_argument('--H', type=float, required=True, help='hauteur')
    parser.add_argument('--alpha', type=float, required=True, help='angle')
    parser.add_argument('--e', type=float, default = 0.05, help='épaisseur')
    parser.add_argument('--R', type=float, default=0.235, help='rayon de la turbine : 0.235')
    parser.add_argument('--L', type=float, default=0.6, help='largeur de la canalisation : 0.6')
    parser.add_argument('--h', type = float, default=0.05, help='meshsize : 0.05')
    parser.add_argument('--simple', action='store_true', help='get a simple deflector if simple')

    parser.add_argument('--output', type=str, default='./output', help='Path to the output directory')

    args = parser.parse_args()

    return args





def get_points(x, H, alpha, e = 0.05, R = 0.235, L = 0.522):
    """
    creates a list of 4 coordinates 
    (those of the points of the rectangular deflector)
    x = deflector position
    H = deflector height
    e = deflector width
    R = turbine radius
    L = pipe diameter
    """
    A = [- x , - L/2]
    B = [- x + e, - L/2]
    C = [- x + e + H/np.abs(np.tan(alpha)), H - L/2]
    D = [- x + H/np.abs(np.tan(alpha)), H - L/2]

    return [A, B, C, D]


def create_simple_deflector(points, h): 
    """
    creates a curve loop for the rectangular deflector
    """


    P, L = [], []
    for p in points : 
        p = gmsh.model.geo.addPoint(p[0], p[1], 0, h)
        P.append(p)
    for i in range(len(P) - 1) :
        l = gmsh.model.geo.addLine(P[i], P[i+1])
        L.append(l)
    L.append(gmsh.model.geo.addLine(P[-1], P[0]))

    return(gmsh.model.geo.addCurveLoops(L)) 


def get_coordinates(C, L):
    """
    gives the end point coordinates 
    of the circle arc for the curved deflector
    """
    R = np.sqrt(C[0]**2 + C[1]**2)
    E_1 = - L/2
    E_0 = - np.sqrt(np.abs(R**2 - E_1**2))
    return [E_0, E_1]


def create_curved_deflector(points, h):
    """
    given a list of points 6 points, 
    and a mesh size h, 
    creates a curveloop for a curved deflector"""
    P = []

    for i in range(len(points)):
        p_i = points[i]
        p = gmsh.model.geo.addPoint(p_i[0], p_i[1], 0, h, i)
        P.append(p)
 

    l1 = gmsh.model.geo.addLine(P[0], P[1])
    l2 = gmsh.model.geo.addLine(P[0], P[3])
    l3 = gmsh.model.geo.addLine(P[2], P[3])
    l4 = gmsh.model.geo.addLine(P[1], P[4])

    circle_arc = gmsh.model.geo.add_circle_arc(P[2], P[5], P[4])
    
    deflector = [l1, l2, l3, l4, circle_arc]


    return(gmsh.model.geo.addCurveLoops(deflector))    


#is the deflector touching/inside the turbine ? 
def is_point_inside_circle(point, center, r):
    distance_squared = (point[0] - center[0])**2 + (point[1] - center[1])**2
    return distance_squared <= r**2

def is_edge_intersect_circle(start_point, end_point, center, r):
    # Check if the line segment intersects the circle
    dx = end_point[0] - start_point[0]
    dy = end_point[1] - start_point[1]

    vec_start_to_center = np.array([center[0] - start_point[0], center[1] - start_point[1]])
    vec_line = np.array([dx, dy])

    # Calculate the projection of the center onto the line
    projection = np.dot(vec_start_to_center, vec_line) / np.dot(vec_line, vec_line)

    # Check if the projection falls within the line segment
    if 0 <= projection <= 1:
        closest_point_on_line = np.array([start_point[0] + projection * dx, start_point[1] + projection * dy])
        distance_squared = np.sum((center - closest_point_on_line)**2)
        return distance_squared <= r**2

    return False

def is_intersection(points_list, r):
    for point in points_list:
        if is_point_inside_circle(point, [0, 0], r):
            return True  # Intersection found

    edges = [(points_list[0], points_list[1]), (points_list[1], points_list[2]),
             (points_list[2], points_list[3]), (points_list[3], points_list[0])]

    for edge in edges:
        if is_edge_intersect_circle(edge[0], edge[1], [0, 0], r):
            return True  # Intersection found

    return False



def generate_filename_prefix(x, H, alpha, simple_defflector):
    # Créer le préfixe du nom de fichier
    if not simple_defflector : 
        filename_prefix = f"curved_deflector_{x:.2f}_{H:.2f}_{alpha:.0f}"
    else : 
        filename_prefix = f"simple_deflector_{x}_{H}_{alpha}"

    # Remplacer les points décimaux par des underscores
    #filename_prefix = filename_prefix.replace('.', '_')

    return filename_prefix


def delete_file(pattern='filename.msh'):
    """
    Supprime tous les fichiers avec le nom spécifié dans le répertoire courant.

    Parameters:
    - pattern (str): Modèle de nom de fichier à supprimer (par défaut 'filename.msh').
    """
    # Utiliser glob pour trouver tous les fichiers correspondant au modèle
    files_to_delete = glob.glob(pattern)

    # Supprimer chaque fichier trouvé
    for file_path in files_to_delete:
        try:
            os.remove(file_path)
            print(f"Fichier supprimé : {file_path}")
        except Exception as e:
            print(f"Erreur lors de la suppression de {file_path}: {e}")


if __name__ == '__main__':

    
    args = get_args()
        # Accéder aux valeurs des arguments
    
    #x = args.x
    # on fixe la position ...
    H = args.H
    alpha = args.alpha
    e = args.e
    R = args.R
    L = args.L
    h = args.h
    simple_deflector = args.simple
    
    output = args.output
    
    """
    # Afficher les valeurs des arguments
    print(f"x: {x}")
    print(f"H: {H}")
    print(f"alpha: {alpha}")
    print(f"e: {e}")
    print(f"R: {R}")
    print(f"L: {L}")
    print(f'h : {h}')
    """

    L, R, e, h = 0.522 , 0.235, 0.05, 0.005
    H = 0.8 * L
    simple_deflector = False

    gmsh.initialize()
    gmsh.model.remove()
    gmsh.model.add("new_model")

    x = R + e + H/np.abs(np.tan(np.radians(alpha)))   
    points = get_points(x, H, np.radians(alpha), e, R, L)

    if is_intersection(points, R) : 
            print('Paramètres non valables : le déflecteur touche la turbine')
    else : 
        if simple_deflector : 
                deflector = create_simple_deflector(points, h)
                gmsh.model.geo.addPlaneSurface(deflector)
        else : 
                C = points[2]
                O = [0, 0]
                E = get_coordinates(C, L)
                points.append(E)
                points.append(O)

                deflector = create_curved_deflector(points, h)
                gmsh.model.geo.addPlaneSurface(deflector)
            
        #filename_prefix = generate_filename_prefix(x, H, alpha, simple_deflector)
        filename_mesh, filename_t = f'{output}.mesh', f'{output}.t'
        
        gmsh.model.geo.synchronize()
        gmsh.model.mesh.generate()
        gmsh.write(filename_mesh)
        gmsh4mtc.gmsh4mtc_single_step(filename_mesh,  filename_t)
        delete_file(filename_mesh) 