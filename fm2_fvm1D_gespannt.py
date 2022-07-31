import numpy as np
import matplotlib.pyplot as plt
import os
import textwrap

##############################
if __name__ == "__main__":
    os.system('clear')

    start_message = "Dieses Programm löst die stationäre Grundwassergleichung in einem gespannten Aquifer in 1D mit der Finite-Volumen-Methode. Vorzugeben sind die Länge des Modellgebiets, die Anzahl der Finiten Volumina (für äquidistante Diskretisierungsweite), die Randbedingungen und mögliche Quell-/Senkterme, sowie die hydraulische Leitfähigkeit."
    print("\n")
    print(textwrap.fill(start_message,80))
    print("\n")

L = input("Länge des Modellgebiets in m: ")
L = float(L)
nCells = input("Anzahl der Zellen (int!): ")
nCells = int(nCells)
Dx = L/nCells
k_f = input("Hydraulische Leitfähigkeit in [m/s]: ")
k_f = float(k_f)
Area = 1 # Querschnittsfläche der Kontrollvolumina einfach mal konstant 1 gesetzt

gitter = np.empty([nCells+1])
gitter[0] = 0.0
for i in range(1,nCells+1):
    gitter[i]= i*Dx

x = gitter
y = np.ones([nCells+1])
plt.figure(dpi=150)
plt.plot(x,y)
plt.grid(axis='x')
plt.yticks([])
plt.xlim([0,L])
plt.ylim([0,1])
plt.xticks(np.arange(min(x), max(x), Dx))
plt.title("Modellgebiet")
plt.show()

coords = np.empty((nCells),float)
node_numbers = np.empty((nCells),int)
coords[0] = 0.0+Dx/2
node_numbers[0]=1
for i in range(1,nCells):
    coords[i] = 0.0+Dx/2+i*Dx
    node_numbers[i] = i+1 # Knotennummer
print("Welche Randbedingungen sollen den Randknoten zugewiesen werden?")
RB_type = np.empty(2,str)
Dirichlet = np.empty(2,float)
Neumann = np.empty(2,float)
RB_type[0] = input("Linker Rand, Dirichlet für Piezometerhöhe (Gib ein: 'D') oder Neumann für Fluss ('N')")
RB_type[1] = input("Rechter Rand, Dirichlet für Piezometerhöhe (Gib ein: 'D') oder Neumann für Fluss ('N')")
if(RB_type[0]=='D'):
    Dirichlet[0] = input("Dirichlet-Piezometerhöhe am linken Rand in m: ")
elif(RB_type[0]=='N'):
    Neumann[0] = input("Neumann-Fluss am linken Rand in m^3/(m^2 s): ")
else:
    print("Typ der Randbedingung am linken Rand nicht korrekt eingegeben. ENDE")
    exit()

if(RB_type[1]=='D'):
    Dirichlet[1] = input("Dirichlet-Piezometerhöhe am rechten Rand in m: ")
elif(RB_type[1]=='N'):
    Neumann[1] = input("Neumann-Fluss am rechten Rand in m^3/(m^2 s): ")
else:
    print("Typ der Randbedingung am rechten Rand nicht korrekt eingegeben. ENDE")
    exit()

A = np.zeros((nCells,nCells),float)
r = np.zeros((nCells),float)
if(RB_type[0]=='D'):
    A[0][0] = 1 # Hauptdiagonale auf 1
    r[0] = Dirichlet[0] # rechte Seite auf den vorgegebenen Dirichlet-Wert
if(RB_type[1]=='D'):
    A[nCells-1][nCells-1] = 1 # Hauptdiagonale auf 1
    r[nCells-1] = Dirichlet[1] # rechte Seite auf den vorgegebenen Dirichlet-Wert
if(RB_type[0]=='N'):
    A[0][0] = -1 # Hauptdiagonale auf -1 (ergibt sich aus der Einführung eines virtuellen Hilfsknotens)
    A[0][1] = 1 # Nachbarknoten bekommt die 1
    r[0] = Neumann[0]*Dx / (-1 * k_f) # berücksichtigen auf der rechten Seite
if(RB_type[1]=='N'):
    A[nCells-1][nCells-1] = -1 # Hauptdiagonale auf -1 (ergibt sich aus der Einführung eines virtuellen Hilfsknotens)
    A[nCells-1][nCells-2] = 1 # Nachbarknoten bekommt die 1
    r[nCells-1] = Neumann[1]*Dx / (-1 * k_f) # berücksichtigen auf der rechten Seite

print("Das Modellgebiet hat die folgenden Knoten (Zellmittelpunkte) mit zugehörigen Koordinaten: ")
for i in range(0,nCells):
    print(node_numbers[i], ":", round(coords[i],2))

# den Quell-/Senkterm einlesen und auf die rechte Seite einbauen
q_pos = input("An welcher Koordinate soll ein Quell-/Senkterm berücksichtigt werden?")
q_pos = float (q_pos)
hilfsarray = coords - q_pos
minwert = abs(hilfsarray).min()
index = 0
for i in range(0,nCells):
    if ((abs(hilfsarray[i]) > minwert-1.e-8) and (abs(hilfsarray[i]) < minwert+1.e-8)):
        index = i
value = input("Zugehöriger Quell-/Senkterm in m^3/s: ")
value = float(value)
r[index] = value*Dx/(-1*k_f*Area)

# die Hauptdiagonale und die Nebeneinträge der Nicht-Randknoten füllen
for i in range(1,nCells-1):
    A[i][i]=-2
    A[i][i-1] = 1
    A[i][i+1] = 1

# jetzt lösen
sol = np.linalg.solve(A,r)

# und dann plotten
x = coords
y = sol
plt.figure(dpi=150)
plt.plot(x,y)
plt.title("Lösung: stationäre Piezometerhöhenverteilung")
plt.xlabel("x [m]")
plt.ylabel("Piezometerhöhe [m]")
plt.show()
