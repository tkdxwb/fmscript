import numpy as np
import matplotlib.pyplot as plt
import os
import textwrap
import pandas as pd
from scipy.interpolate import griddata
from matplotlib import cm
from matplotlib.ticker import LinearLocator,FormatStrFormatter

### Klasse für die Knoten
class Knoten:
    def __init__(self,nummer,volumen,x_coord,y_coord,n_nachbar,s_nachbar,w_nachbar,o_nachbar,flag_rand_intern):
        self.nummer = nummer
        self.volumen = volumen
        self.x_coord = x_coord
        self.y_coord = y_coord
        self.n_nachbar = n_nachbar
        self.s_nachbar = s_nachbar
        self.w_nachbar = w_nachbar
        self.o_nachbar = o_nachbar
        self.flag = flag_rand_intern

    def __print_Knoten_info_alle__(self):
        print("Nummer:", self.nummer,"x,y:", self.x_coord,self.y_coord, "Nachbarknoten n,s,w,o:", self.n_nachbar,self.s_nachbar,self.w_nachbar,self.o_nachbar)

    def __print_Knoten_info_rand__(self):
        if(self.flag==0):
            print("Nummer:", self.nummer,"x,y:", self.x_coord,self.y_coord, "Nachbarknoten n,s,w,o:", self.n_nachbar,self.s_nachbar,self.w_nachbar,self.o_nachbar)

    def __print_Knoten_info_intern__(self):
        if(self.flag==1):
            print("Nummer:", self.nummer,"x,y:", self.x_coord,self.y_coord, "Nachbarknoten n,s,w,o:", self.n_nachbar,self.s_nachbar,self.w_nachbar,self.o_nachbar)


def berechnen(L, B, nxCells, nyCells, k_f, depth, isJupyter = True):
    Dx = L/nxCells
    Dy = B/nyCells

    gitterx = np.empty([nxCells+1])
    gitterx[0] = 0.0
    for i in range(1,nxCells+1):
        gitterx[i]= i*Dx
    gittery = np.empty([nxCells+1])
    gittery.fill(B)

    x = gitterx
    y = gittery

    plt.ioff()
    figure = plt.figure(dpi=150)
    plt.plot(x,y)
    plt.grid()
    plt.xlim([0,L])
    plt.ylim([0,B])
    plt.xticks(np.arange(min(x), max(x), Dx),color='w')
    plt.yticks(np.arange(0, max(y), Dy),color='w')
    plt.title("Modellgebiet zur Kontrolle, Fenster schließen um weiter zu machen")
    display(figure) if isJupyter else plt.show()


    number_of_cells = nxCells*nyCells
    list = []
    volumen = Dx*Dy*depth
    rowcounter = 1
    for i in range(1,number_of_cells+1):
        ### unterste Reihe
        rowcounter = (i-1)//nxCells + 1
        #print("i, rowcounter:",i,rowcounter)
        if (i==1):
            x = Dx/2
            y = Dy/2
            flag = 0
            n = nxCells+1
            o = i+1
            w = 0
            s = 0
            #print(n,s,w,o)
        if ((i>1) and (i<nxCells)):
            x = (i-1)*Dx + Dx/2
            y = Dy/2
            flag = 0
            n = nxCells+i
            o = i+1
            w = i-1
            s = 0
            #print(n,s,w,o)
        if (i==nxCells):
            x = (i-1)*Dx + Dx/2
            y = Dy/2
            flag = 0
            n = nxCells+i
            o = 0
            s = 0
            w = i-1
            #print(n,s,w,o)
        ### die Reihen dazwischen ...
        if((i>nxCells) and (i==((rowcounter-1)*nxCells+1))):
            x = Dx/2
            y = (rowcounter-1)*Dy + Dy/2
            flag = 0
            n = i+nxCells
            o = i+1
            s = i-nxCells
            w = 0
        if((i>nxCells) and (i>((rowcounter-1)*nxCells+1)) and (i<nxCells*rowcounter)):
            x = (i-1-(rowcounter-1)*nxCells)*Dx + Dx/2
            y = (rowcounter-1)*Dy + Dy/2
            flag = 1 # das sind die internen Knoten!
            n = i+nxCells
            o = i+1
            s = i-nxCells
            w = i-1
        if(i==nxCells*rowcounter):
            x = (i-1-(rowcounter-1)*nxCells)*Dx + Dx/2
            y = (rowcounter-1)*Dy + Dy/2
            flag = 0
            n = i+nxCells
            o = 0
            s = i-nxCells
            w = i-1
        ### oberste Reihe
        if(i==(number_of_cells-nxCells+1)):
            x = Dx/2
            y = (rowcounter-1)*Dy + Dy/2
            flag = 0
            n = 0
            o = i+1
            s = i-nxCells
            w = 0
        if((i>(number_of_cells-nxCells+1)) and (i<number_of_cells)):
            x = (i-1-(rowcounter-1)*nxCells)*Dx + Dx/2
            y = (rowcounter-1)*Dy + Dy/2
            flag = 0
            n = 0
            o = i+1
            s = i-nxCells
            w = i-1
        if(i==number_of_cells):
            x = (i-1-(rowcounter-1)*nxCells)*Dx + Dx/2
            y = (rowcounter-1)*Dy + Dy/2
            flag = 0
            n = 0
            o = 0
            s = i-nxCells
            w = i-1
        ### Knoten der List hinzufügen
        list.append(Knoten(i,volumen,x,y,n,s,w,o,flag))

    ### Kontrollausgaben, kann wieder auskommentiert werden
    #for node in list:
        #node.__print_Knoten_info_alle__()
    #
    #print("\n", "Randknoten")
    #
    #for node in list:
        #node.__print_Knoten_info_rand__()
    #
    #print("\n", "Interne Knoten")
    #
    #for node in list:
        #node.__print_Knoten_info_intern__()

    ### Randbedingungen einlesen
    message = "Jetzt müssen Randbedingungen zugewiesen werden. Im Moment erlaubt die Implementierung eine Zuweisung nur pro Randabschnitt (nord, süd, west, ost), entweder Dirichlet oder Neumann."
    print("\n")
    print(textwrap.fill(message,80))
    print("\n")

    RB_type = np.empty(4,str)
    Dirichlet = np.empty(4,float)
    Neumann = np.empty(4,float)
    RB_type[0] = input("Linker Rand (west), Dirichlet für Piezometerhöhe (Gib ein: 'D') oder Neumann für Fluss ('N')")
    RB_type[1] = input("Rechter Rand (ost), Dirichlet für Piezometerhöhe (Gib ein: 'D') oder Neumann für Fluss ('N')")
    RB_type[2] = input("Oberer Rand (nord), Dirichlet für Piezometerhöhe (Gib ein: 'D') oder Neumann für Fluss ('N')")
    RB_type[3] = input("Unterer Rand (süd), Dirichlet für Piezometerhöhe (Gib ein: 'D') oder Neumann für Fluss ('N')")
    if(RB_type[0]=='D'):
        Dirichlet[0] = input("Dirichlet-Piezometerhöhe am linken Rand (west) in m: ")
    elif(RB_type[0]=='N'):
        Neumann[0] = input("Neumann-Fluss am linken Rand (west) in m^3/(m^2 s): ")
    else:
        print("Typ der Randbedingung am linken Rand (west) nicht korrekt eingegeben. ENDE")
        exit()
    if(RB_type[1]=='D'):
        Dirichlet[1] = input("Dirichlet-Piezometerhöhe am rechten Rand (ost) in m: ")
    elif(RB_type[1]=='N'):
        Neumann[1] = input("Neumann-Fluss am rechten Rand (ost) in m^3/(m^2 s): ")
    else:
        print("Typ der Randbedingung am rechten Rand (ost) nicht korrekt eingegeben. ENDE")
        exit()
    if(RB_type[2]=='D'):
        Dirichlet[2] = input("Dirichlet-Piezometerhöhe am oberen Rand (nord) in m: ")
    elif(RB_type[2]=='N'):
        Neumann[2] = input("Neumann-Fluss am oberen Rand (nord) in m^3/(m^2 s): ")
    else:
        print("Typ der Randbedingung am oberen Rand (ost) nicht korrekt eingegeben. ENDE")
        exit()
    if(RB_type[3]=='D'):
        Dirichlet[3] = input("Dirichlet-Piezometerhöhe am unteren Rand (süd) in m: ")
    elif(RB_type[3]=='N'):
        Neumann[3] = input("Neumann-Fluss am unteren Rand (süd) in m^3/(m^2 s): ")
    else:
        print("Typ der Randbedingung am unteren Rand (süd) nicht korrekt eingegeben. ENDE")
        exit()

    ### Quell-/Senkterm einlesen
    q_posx = input("An welcher x-Koordinate soll ein Quell-/Senkterm berücksichtigt werden?")
    q_posx = float (q_posx)
    q_posy = input("An welcher y-Koordinate soll ein Quell-/Senkterm berücksichtigt werden?")
    q_posy = float (q_posy)
    q_term = input("Quell-/Senkterm in m^3/s):")
    q_term = float (q_term)
    Aquifer_Tiefe = 1 ### einfach mal gleich eins gesetzt
    hilfsarray = np.empty((0,2),dtype=np.float64)
    i= 0
    for node in list:
        dx = abs(node.x_coord - q_posx)
        dy = abs(node.y_coord - q_posy)
        neue_Zeile = np.array([[dx,dy]], dtype=np.float64)
        hilfsarray = np.append(hilfsarray,neue_Zeile,axis=0)
        i = i+1
    no_of_nodes = i

    minInColumns = np.amin(hilfsarray, axis=0)
    #print(hilfsarray)
    #print(minInColumns)
    #print(no_of_nodes)
    index = 0
    for node in list:
        if ((abs(node.x_coord-q_posx)>minInColumns[0]-1.e-4) and (abs(node.x_coord-q_posx)<minInColumns[1]+1.e-4) and (abs(node.y_coord-q_posy)>minInColumns[1]-1.e-4) and (abs(node.y_coord-q_posy)<minInColumns[1]+1.e-4)):
            index = node.nummer
    #print("Knotennummer für Quellterm ist: ", index)

    ### Jetzt beginnen mit dem Zusammenbauen der Matrix und der rechten Seite

    A = np.zeros((no_of_nodes,no_of_nodes),float)
    r = np.zeros((no_of_nodes),float)
    ### Zunächst die Randbedingungen
    for node in list:
        if ((node.x_coord<Dx) and (node.y_coord>Dy) and (node.y_coord<(B-Dy))): ### dann sind wir auf dem linken Rand (Eckknoten gehen extra)
            if (RB_type[0]=='D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[0]
            if (RB_type[0]=='N'):
                A[node.nummer-1][node.nummer-1] = 2*(Dx/Dy+Dy/Dx) - Dy/Dx
                A[node.nummer-1][node.o_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.s_nachbar-1] = -1*Dx/Dy
                A[node.nummer-1][node.n_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*Neumann[0]/k_f
        if ((node.x_coord>(L-Dx)) and (node.y_coord>Dy) and (node.y_coord<(B-Dy))): ### dann sind wir auf dem rechten Rand (Eckknoten gehen extra)
            if (RB_type[1]=='D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[1]
            if (RB_type[1]=='N'):
                A[node.nummer-1][node.nummer-1] = 2*(Dx/Dy+Dy/Dx) - Dy/Dx
                A[node.nummer-1][node.w_nachbar-1] = -1*Dy/Dx 
                A[node.nummer-1][node.s_nachbar-1] = -1*Dx/Dy
                A[node.nummer-1][node.n_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*Neumann[1]/k_f
        if ((node.y_coord<Dy) and (node.x_coord>Dx) and (node.x_coord<(L-Dx))): ### unterer Rand (Eckknoten gehen extra)
            if (RB_type[3]=='D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[3]
            if (RB_type[3]=='N'):
                A[node.nummer-1][node.nummer-1] = 2*(Dx/Dy+Dy/Dx) - Dx/Dy
                A[node.nummer-1][node.o_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.w_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.n_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*Neumann[3]/k_f
        if ((node.y_coord>(B-Dy)) and (node.x_coord>Dx) and (node.x_coord<(L-Dx))): ### oberer Rand (Eckknoten gehen extra)
            if (RB_type[2]=='D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[2]
            if (RB_type[2]=='N'):
                A[node.nummer-1][node.nummer-1] = 2*(Dx/Dy+Dy/Dx) - Dx/Dy
                A[node.nummer-1][node.o_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.w_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.s_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*Neumann[2]/k_f
        if ((node.y_coord<Dy) and (node.x_coord<Dx)): ### Eckknoten unten links
            if (RB_type[0]=='D'): 
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[0]
            if (RB_type[3] == 'D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[3]
            if ((RB_type[0]=='N') and (RB_type[3] == 'N')): # wenn beide Neumann sind muss der Eckknoten auch Neumann werden
                print("links unten beide Neumann")
                A[node.nummer-1][node.nummer-1] = Dx/Dy+Dy/Dx
                A[node.nummer-1][node.o_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.n_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*(Neumann[3]+Neumann[0])/k_f
        if ((node.y_coord<Dy) and (node.x_coord>(L-Dx))): ### Eckknoten unten rechts
            if (RB_type[1]=='D'): 
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[1]
            if (RB_type[3] == 'D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[3]
            if ((RB_type[1]=='N') and (RB_type[3] == 'N')): # wenn beide Neumann sind muss der Eckknoten auch Neumann werden
                print("rechts unten beide Neumann")
                A[node.nummer-1][node.nummer-1] = Dx/Dy+Dy/Dx
                A[node.nummer-1][node.w_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.n_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*(Neumann[3]+Neumann[1])/k_f
        if ((node.y_coord>(B-Dy)) and (node.x_coord<Dx)): ### Eckknoten oben links
            if (RB_type[0]=='D'): 
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[0]
            if (RB_type[2] == 'D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[2]
            if ((RB_type[0]=='N') and (RB_type[2] == 'N')): # wenn beide Neumann sind muss der Eckknoten auch Neumann werden
                print("links oben beide Neumann")
                A[node.nummer-1][node.nummer-1] = Dx/Dy+Dy/Dx
                A[node.nummer-1][node.o_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.s_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*(Neumann[2]+Neumann[0])/k_f
        if ((node.y_coord>(B-Dy)) and (node.x_coord>(L-Dx))): ### Eckknoten oben rechts
            if (RB_type[1]=='D'): 
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[1]
            if (RB_type[2] == 'D'):
                A[node.nummer-1][node.nummer-1] = 1
                r[node.nummer-1] = Dirichlet[2]
            if ((RB_type[1]=='N') and (RB_type[2] == 'N')): # wenn beide Neumann sind muss der Eckknoten auch Neumann werden
                print("rechts oben beide Neumann")
                A[node.nummer-1][node.nummer-1] = Dx/Dy+Dy/Dx
                A[node.nummer-1][node.w_nachbar-1] = -1*Dy/Dx
                A[node.nummer-1][node.s_nachbar-1] = -1*Dx/Dy
                r[node.nummer-1] = -1*(Neumann[2]+Neumann[1])/k_f
        if (node.flag==1): ### alle internen Knoten sind einfach
            A[node.nummer-1][node.nummer-1] = 2*(Dx/Dy+Dy/Dy)
            A[node.nummer-1][node.n_nachbar-1] = -1*Dx/Dy
            A[node.nummer-1][node.s_nachbar-1] = -1*Dx/Dy
            A[node.nummer-1][node.w_nachbar-1] = -1*Dy/Dx
            A[node.nummer-1][node.o_nachbar-1] = -1*Dy/Dx
            r[node.nummer-1] = 0
        if (node.nummer==index): ### Quell-/Senkterm dazu
            r[node.nummer-1] = r[node.nummer-1] - q_term/(k_f*Aquifer_Tiefe)

    #print("Matrix:")
    #print(A)
    #print("rechte Seite:")
    #print(r)

    ### jetzt lösen
    sol = np.linalg.solve(A,r)
    #print("Lösung:")
    #print(sol)

    ### jetzt plotten
    data = np.zeros((no_of_nodes,3),float)
    i = 0
    for node in list:
        data[node.nummer-1][0] = node.x_coord
        data[node.nummer-1][1] = node.y_coord
        data[node.nummer-1][2] = sol[node.nummer-1]

    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    xyz = {'x': x, 'y': y, 'z': z}

    # put the data into a pandas DataFrame
    df = pd.DataFrame(xyz, index=range(len(xyz['x'])))

    # re-create the 2D-arrays
    x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
    y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
    x2, y2 = np.meshgrid(x1, y1)
    z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')

    fig = plt.figure(dpi=150)
    ax = fig.add_subplot(projection='3d')
    surf = ax.plot_surface(x2, y2, z2, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    #ax.set_zlim(0, 100)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title('Plot der Piezometerhöhen')

    fig.show()


##############################
if __name__ == "__main__":
    os.system('clear')

    start_message = "Dieses Programm löst die stationäre Grundwassergleichung in einem gespannten Aquifer in 2D mit der Finite-Volumen-Methode. Vorzugeben sind die Länge des Modellgebiets, die Anzahl der Finiten Volumina in x- und y-Richtung (für jeweils äquidistante Diskretisierungsweite), die Randbedingungen und mögliche Quell-/Senkterme, sowie die hydraulische Leitfähigkeit."
    print("\n")
    print(textwrap.fill(start_message,80))
    print("\n")

    L = input("Länge des Modellgebiets in m: ")
    L = float(L)
    B = input("Breite des Modellgebiets in m: ")
    B = float(B)
    nxCells = input("Anzahl der Zellen in x-Richtung (int!): ")
    nxCells = int(nxCells)
    nyCells = input("Anzahl der Zellen in y-Richtung (int!): ")
    nyCells = int(nyCells)
    Dx = L/nxCells
    Dy = B/nyCells
    k_f = input("Hydraulische Leitfähigkeit in [m/s]: ")
    k_f = float(k_f)
    depth = 1 # Querschnittstiefe der Kontrollvolumina einfach mal konstant 1 gesetzt

    print("Die Nummerierung der Finiten Volumina erfolgt, beginnend mit '1' bei (0/0) zunächst in positive x-Richtung, dann in positive y-Richtung.")

    berechnen(L,B,nxCells, nyCells, k_f, depth, isJupyter = False)