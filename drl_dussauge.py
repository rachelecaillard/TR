#--- GENERIC IMPORT ---------------------------------------------------------------------------------------------------------+
import os
import sys
import math
import time
import numpy as np
import gmsh
import matplotlib as plt
import datetime as dt
import create_deflector as cd
import analyse_and_reward as aar
import gmsh4mtc



        #################################################
        ## ENVIRONMENT TO OPTIMISE AN AIRFOIL WITH DRL ##
        #################################################


class drl_dussauge():

#--- CREATE OBJECT ---------------------------------------------------------------------------------------------------------+
    def __init__(self, path):
        self.name     = 'drl_dussauge'                                     
        self.act_size = 2
        self.obs_size = self.act_size
        self.obs      = np.zeros(self.obs_size)

        # attention alpha en degrés
        self.x_min    = np.array([0.261, 10])  # (H, alpha) min 
        self.x_max    = np.array([0.46, 85])   # (H, alpha) max
        self.x_0      = 0.5*(self.x_min+self.x_max)

        self.path        = path                                                       
        self.episode     = 0                   # Set episode number

        self.time_init   = 0.                  # Pour connaître la durée d'un épisode
        self.time_end    = 0.


    #--- METHODES DE LA CLASSE ----------------------------------------------------------------------------------------------+

    def solve_problem_cimlib(self):
        """ Solve problem using cimlib and move vtu and drag folder. It changes properties."""

        os.system('cd '+self.output_path+'cfd/.; touch run.lock; mpirun -n 8 /softs/cemef/cimlibxx/master/bin/cimlib_CFD_driver'+' Principale.mtc > trash.txt;')
        
        #on déplace les vtu
        os.system('mv '+self.output_path+'cfd/Resultats/2d/* '+self.vtu_path+'.')

        #on déplace les fichiers dans le dossier data
        os.system( 'mv ' + self.output_path + 'cfd/Resultats/Torque.txt ' + self.data + '.')
        os.system('mv '+self.output_path+'cfd/Resultats/capteurP1.txt '+ self.data + '.')
        os.system('mv '+self.output_path+'cfd/Resultats/capteurP2.txt '+ self.data + '.')
        os.system('mv '+self.output_path+'cfd/Resultats/Variables.txt '+ self.data + '.')

        os.system('rm -r '+self.output_path+'cfd')
        os.system('cp -r '+self.vtu_path+'bulles_00500.vtu ./video/')                                  
        os.system('mv ./video/bulles_00500.vtu '+'./video/video_'+str(self.episode)+'.vtu')
        

    def shape_generation_deflector(self, x, mesh_path, t_path):
        """ Generate shape """
        H, alpha = x[0], x[1]
        mesh_size      = 0.005                                                
        try:
            L, R, e = 0.522 , 0.235, 0.05
            H = 0.8 * L
            simple_deflector = False

            gmsh.initialize()
            """gmsh.model.remove()
            gmsh.model.add("new_model")"""

            x = R + e + H/np.abs(np.tan(np.radians(alpha)))   
            points = cd.get_points(x, H, np.radians(alpha), e, R, L)

            if cd.is_intersection(points, R) : 
                    print('Paramètres non valables : le déflecteur touche la turbine')
            else : 
                if simple_deflector : 
                        deflector = cd.create_simple_deflector(points, mesh_size)
                        gmsh.model.geo.addPlaneSurface(deflector)
                else : 
                        C = points[2]
                        O = [0, 0]
                        E = cd.get_coordinates(C, L)
                        points.append(E)
                        points.append(O)

                        deflector = cd.create_curved_deflector(points, mesh_size)
                        gmsh.model.geo.addPlaneSurface(deflector)
                
                gmsh.model.geo.synchronize()
                gmsh.model.mesh.generate()
                gmsh.write(mesh_path)
                
                gmsh4mtc.gmsh4mtc_single_step(mesh_path, t_path)
                cd.delete_file(mesh_path) 
                
                

        except Exception as e:
            pass
           

    def compute_reward(self, x, path_folder):
        """ Calcule le reward puis l'enregistre dans un Values.txt avec les autres paramètres"""
        Tau, Tau_f = 12, 18
        deltaP_0 = 10000
        reward = aar.get_reward(path_folder, deltaP_0, Tau, Tau_f)
        if reward is not None :
            self.reward = reward['Reward'].mean()
        else :
            self.reward = -1000


        ### Ecriture dans Values
        """print(os.path)
        if not os.path.isfile('Values.txt'):
            f = open('Values.txt','w')
            f.write(
                'Index'+'\t'+'H'+'\t'+'alpha'+'\t'+'Reward'+'\n'
                )
        else:
            f = open('Values.txt','a')
        f.write(
                str(self.episode)+'\t'+"{:.3e}".format(x[0])+ 't'+"{:.3e}".format(x[1])+'\t'"{:.3e}".format(self.reward)+'\n'
            )
        f.close()"""
        
        self.episode += 1 


    #--- CFD RESOLUTION -------------------------------------------------------------------------------------------------+

    def cfd_solve(self, x, ep):
        """
        lance la simulation, 
        crée la forme du déflecteur, 
        convertit le maillage en .t, 
        résout la simulation, 
        enregistre la récompense
        """

        self.time_init=dt.datetime.now()                                        # On suit en temps le DRL
        if not os.path.isfile('temps_start.txt'):
            f = open('temps_start.txt','w')
            f.write('Index'+'\t'+'Heure start'+'\n')
            f.close()
        f = open('temps_start.txt','a')
        f.write(str(ep)+'\t'+ dt.datetime.now().strftime("%H:%M:%S")+'\n')
        f.close()

        ### Create folders and copy cfd (please kill me)
        ### On met les résultats là dedans 
        self.output_path = self.path+'/'+str(ep)+'/'  # Pour chaque épisode
        self.vtu_path    = self.output_path+'vtu/'
        self.data = self.output_path+'data/'
        self.msh_path    = self.output_path+'msh/'
        self.t_mesh_path = self.output_path+'t_mesh/'
        
        os.makedirs(self.data)
        os.makedirs(self.vtu_path)
        os.makedirs(self.msh_path)
        os.makedirs(self.t_mesh_path)
        os.system('cp -r cfd ' + self.output_path + '.')   


        ### create the shape 
        self.shape_generation_deflector(x, self.output_path + 'cfd/deflector.msh', self.output_path + 'cfd/deflector.t')                        

        """ ### convert .mesh to .t
        os.system('cd '+self.output_path+'cfd ; python3 gmsh2mtc.py')
        os.system('cd '+self.output_path+'cfd ; cp -r deflector.msh ../msh')
        os.system('cd '+self.output_path+'cfd ; module load cimlibxx/master')
        os.system('cd '+self.output_path+'cfd ; echo 0 | mtcexe deflector.t')
        os.system('cd '+self.output_path+'cfd ; cp -r deflector.t ../t_mesh')"""
        
        ### solving the problem
        self.solve_problem_cimlib()

        ### Compute the reward 
        #self.compute_reward(self.data)
        self.reward = 0

        ### On écrit la durée
        self.time_end     = dt.datetime.now()
        difference        = self.time_end - self.time_init
        heures, reste     = divmod(difference.seconds, 3600)
        minutes, secondes = divmod(reste, 60)
        
        """if not os.path.isfile('duree.txt'):
            f = open('duree.txt','w')
            f.write('Index'+'\t'+'Heure start'+'\t'+'Heure end'+'\t'+'Durée'+'\n')
            f.close()
        fi = open('duree.txt','a')
        fi.write(
            str(ep)+'\t'+ self.time_init.strftime("%H:%M:%S")+'\t'
            +self.time_end.strftime("%H:%M:%S")+'\t'
            +f"{str(heures)}:{str(minutes)}:{str(secondes)}"+'\n'
            )
        fi.close()"""
        return self.reward

                

    ### Take one step
    def step(self, actions, ep):
        conv_actions = self.convert_actions(actions)
        reward       = self.cfd_solve(conv_actions, ep)
        return reward, conv_actions

    ### Provide observation
    def observe(self):
        # Always return the same observation
        return self.obs

    ### Convert actions
    def convert_actions(self, actions):
        """ Convertit les actions du DRL qui sont entre 0 et 1 """
        # Convert actions
        conv_actions  = self.act_size*[None]
        x_p           = self.x_max - self.x_0
        x_m           = self.x_0   - self.x_min

        for i in range(self.act_size):
            if (actions[i] >= 0.0):
                conv_actions[i] = self.x_0[i] + x_p[i]*actions[i]
            if (actions[i] <  0.0):
                conv_actions[i] = self.x_0[i] + x_m[i]*actions[i]
        return conv_actions

    ### Close environment
    def close(self):
        pass

    ### A function to replace text in files
    ### This function finds line containing string, erases the
    ### whole line it and replaces it with target
    def line_replace(self, string, line, target):
        command = "sed -i '/"+string+"/c\\"+line+"' "+target
        os.system(command)


#--- FONCTIONS DE PARAMETRISATION ----------------------------------------------------------------------------------------------+
   

#--- FONCTION DE PUNITION POUR LA SURFACE -------------------------------------------------------------------------------------+

    def punition_exponentielle(self):
        """ Donne la punition que l'on doit mettre dans le reward (exponentielle) """
        if self.area < self.area_min :
            return np.exp((self.area_min/self.area) -1) - 1       # vaut 0 au début 
        else : 
            return 0. 

    def punition_affine_marge(self, marge):
        """ Donne une punition affine de alpha * (S-Sref) avec marge de marge %"""
        if self.area_target * (1 - marge) < self.area < self.area_target * (1 + marge) : 
            return 0
        elif self.area < self.area_target * (1 - marge) : 
            return self.alpha * abs((1 - marge) * self.area_target - self.area)
        else : 
            return self.alpha * abs((1 + marge) * self.area_target - self.area)

    def punition_affine(self) : 
        """ Donne une punition de la forme alpha * abs(S-Sref) """
        return self.alpha * abs(self.area - self.area_target)
