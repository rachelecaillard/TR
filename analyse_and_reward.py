import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os



def read_data(path_folder):
    path_torque = os.path.join(path_folder, 'Torque.txt')
    path_variables = os.path.join(path_folder, 'Variables.txt')
    path_capteurP1 = os.path.join(path_folder, 'capteurP1.txt')
    path_capteurP2 = os.path.join(path_folder, 'capteurP2.txt')

    df_t = pd.read_csv(path_torque, sep="\t", header=0, usecols=[0, 1, 2])
    param = pd.read_csv(path_variables, sep="\t", header=0, usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8])
    df_P1 = pd.read_csv(path_capteurP1, sep="\t", header=0, usecols=[0, 1])
    df_P2 = pd.read_csv(path_capteurP2, sep="\t", header=0, usecols=[0, 1])
    df_P = df_P1 - df_P2
    df_P['Temps'] = df_P1['Temps']
    return df_t, param, df_P


def show_param(param, C_ech = 1):
    R = param.loc[0]['R_turb']/C_ech
    Rho = 1000 #param.loc[0]['Rho1']
    Eta = param.loc[0]['Eta1']
    Veau = param.loc[0]['VitesseAir']
    Vrot = param.loc[0]['VitesseRotation']
    print("Rayon turbine =", R, 'm')
    print("Rho =", Rho, 'kg/m^3')
    print("Eta =", Eta, 'Pa.s')
    print("Vitesse air =", Veau, 'm/s')
    print("Vitesse rotation = ", Vrot, 'rad/s')

    Re = 2*Rho*Veau*R/Eta
    print("Reynolds = ", Re)

    Lambda = R*Vrot/Veau
    print("Tip speed ratio = ", Lambda)


def delete_transi(df, Tau = 1, Tau_f = 0): 

    iNonValide = df[df['Temps']<Tau].index.tolist()
    iDepart = iNonValide[-1]

    if Tau_f < len(df) : 
        iNonValide_2 =  df[df['Temps'] > Tau_f].index.tolist()
        iFin = iNonValide_2[0]
        
        if not iNonValide_2:
            iFin = len(df)
        else:
            iFin = iNonValide_2[0]
        
        return df[iDepart:iFin]

    return df[iDepart :]


def get_Cp(df, Vrot, Rho, Veau) : 
    Cp = pd.DataFrame()
    Cp['Cp'] = df['Torque']*Vrot/(0.5*Rho*Veau**3)
    Cp['Temps'] = df['Temps']
    return Cp


def plot_analysis(Cp, df_P, r) : 
    fig, axs = plt.subplots(3, 1, figsize=(5, 4))

    Cp_moy = Cp["Cp"].mean()

    axs[0].plot(Cp['Temps'], Cp['Cp'])
    axs[0].set_title('Cp(t), ' f'$\\bar{{Cp}} = {Cp_moy:.3f}$')

    P_moy = df_P['Pression'].mean()

    axs[1].plot(df_P['Temps'], df_P['Pression'])
    axs[1].set_title('Perte de charge (t), ' f'$\\bar{{\Delta P}} = {P_moy:.3f}$')

    axs[2].plot(r['Temps'], r['Reward'])
    axs[2].set_title('récompense(t)')


    plt.tight_layout()

    #plt.title('$\bar{r} = $' + str(r['Reward'].mean()))
    plt.title('r = Cp - $\Delta P/10^4 $, 'f'$\\bar{{r}} = {r["Reward"].mean():.3f}$')

    plt.show()


def get_reward(path_folder,  deltaP_0, Tau = 2, Tau_f = 0,C_ech = 1):
    df_t, param, df_P = read_data(path_folder)
    if Tau_f == None or Tau_f > len(df_t):
        Tau_f = len(df_t)
    
    R, Rho, Eta, Veau, Vrot = param.loc[0]['R_turb']/C_ech, param.loc[0]['Rho1'], param.loc[0]['Eta1'], param.loc[0]['VitesseAir'], param.loc[0]['VitesseRotation']
    Rho = 1000 #je ne sais pas où modifier ça dans les mtc ...
    R = 0.235
    show_param(param, C_ech)
    df_t, df_P = delete_transi(df_t, Tau, Tau_f), delete_transi(df_P, Tau, Tau_f)

    Cp = get_Cp(df_t, Vrot, Rho, Veau)  

    r = pd.DataFrame()
    r['Temps'] = df_P['Temps']
    r['Reward'] = Cp['Cp'] - (df_P['Pression'] - deltaP_0)/deltaP_0

    print('Cp =',  Cp['Cp'].mean())
    print('Delta P =', df_P['Pression'].mean())

    #plot_analysis(Cp, df_P, r)

    return r
    
    
if __name__ == "__main__":
    path_folder = sys.argv[1]
    Tau = 12
    Tau_f = 18
    deltaP_0 = 10000
    r = get_reward(path_folder, Tau, Tau_f)
    print(r['Reward'].mean())