from cProfile import label
from colorsys import hls_to_rgb
import numpy as np
import math
import matplotlib.pyplot as plt
import textwrap
import os

### Beginn: Zur Berechnung der DGL benötigte Funktionen
def dh_von_dx(dx,I_0,Q,k_str,r_hy,A,h):
    zaehler = dx * (I_0 - Q*Q/(k_str * k_str * r_hy**(4/3) * A*A))
    nenner = 1. - Q*Q/(9.81*h*A*A)
    dh = zaehler/nenner
    return dh

def berechne_querschnitt(h,b,m):
    A = b*h + m*h*h
    L_u = b + 2*h*math.sqrt(1.+m*m)
    r_hy = A/L_u
    return(A,r_hy)
### Ende: Zur Berechnung der DGL benötigte Funktionen


### Beginn: Berechnung des Normalabflusszustands und Grenzabfluss
def berechne_f_y (y0,b,I0,Q,kStr,m):
    A, r_hy = berechne_querschnitt(y0,b,m)
    f_y = (A * kStr * r_hy**(2/3) * I0**0.5) - Q
    return f_y

def berechne_normalabfluss(Q,b,I_0,k_str,m):
    y0 = Q/b # Startwertschätzung für Iteration (gibt's bessere Ideen?)
    yneu = y0
    yalt = y0+1 # um die Schleife sicher zu aktivieren
    while abs(berechne_f_y(yneu,b,I_0,Q,k_str,m))>0.001:
        f_y = berechne_f_y (yneu,b,I_0,Q,k_str,m)
        # Ableitung numerisch berechnen
        dy = 0.00000001*f_y
        fplus = berechne_f_y (yneu+dy,b,I_0,Q,k_str,m)
        fminus = berechne_f_y (yneu-dy,b,I_0,Q,k_str,m)
        fabl = (fplus-fminus)/(2*dy)
        yalt = yneu
        #print("aktueller Iterationswert:", yneu, "\n")
        yneu = yneu - (1/fabl)*f_y
    abs(berechne_f_y(yneu,b,I_0,Q,k_str,m))
    h_N = yneu
    v_N = Q/(yneu*b + m*yneu*yneu)
    Fr_N = v_N/(9.81*yneu)**0.5
    # jetzt h_gr ausrechnen (bissle murksige Methode in zwei Schritten weil Trapez)
    q = Q/b
    h_gr = (q*q/9.81)**(1/3)
    q = Q/(b+m*h_gr)
    h_gr = (q*q/9.81)**(1/3)

    return (h_N,v_N,Fr_N,h_gr)
### Ende: Berechnung des Normalabflusszustands und Grenzabfluss

### Beginn: Definition der Klasse Spiegellinie
class Spiegellinie:
    def __init__(self, Q,b,I_0,k_str,m, name):
        self.Q = Q
        self.b = b
        self.I_0 = I_0
        self.k_str = k_str
        self.m = m
        self.h_N, self.v_N, self.Fr_N, self.h_gr = berechne_normalabfluss(Q,b,I_0,k_str,m)
        self.name = name

    def __startWerte__(self, h_start, grenze):
        self.h_start = h_start
        self.grenze = grenze

    def __dxWerte__(self, dx, dx_scale):
        self.dx = dx
        self.dx_scale = dx_scale

    def normalAbfluss(self):
        print("Berechnet wurden: ")
        print("Normalabflusstiefe:", round(self.h_N,2), "m")
        print("Geschwindigkeit bei Normalabfluss:", round(self.v_N,2), "m/s")
        print("Froude-Zahl:", round(self.Fr_N,2))
        print("H_0 bei Normalabfluss:", round(self.h_N+self.v_N**2/19.62,2), "m")
        print("Grenzabflusstiefe:", round(self.h_gr,2), "m")

    def berechnung(self):
        print("Berechnung der Spiegellinie.")
        h_0 = self.h_start
        z_0 = 0.0
        z_So = z_0
        z_Wsp = z_So+h_0
        A,r_hy = berechne_querschnitt(h_0,self.b,self.m)
        v = self.Q/A
        z_H = z_So + h_0 + v*v/19.62
        x = 0
        dx = self.dx
        kurven = np.empty((0,4),dtype=np.float64)
        neue_Zeile = np.array([[x,z_So,z_Wsp,z_H]], dtype=np.float64)
        kurven = np.append(kurven,neue_Zeile,axis=0)
        #print(kurven)
        while (abs(h_0-self.grenze)>0.01):
            x = x+dx
            A,r_hy = berechne_querschnitt(h_0,self.b,self.m)
            v = self.Q/A
            dh = dh_von_dx(dx,self.I_0,self.Q,self.k_str,r_hy,A,h_0)
            h_0 = h_0 + dh
            z_So = z_0 - x*self.I_0
            z_Wsp = z_So + h_0
            z_H = z_So + h_0 + v*v/19.62
            neue_Zeile = np.array([[x,z_So,z_Wsp,z_H]], dtype=np.float64)
            kurven = np.append(kurven,neue_Zeile,axis=0)
            #print(kurven)
            dx = dx * self.dx_scale
            # jetzt die Plots
        x_plot = kurven[:,0]
        y_plot_1 = kurven[:,1]
        y_plot_2 = kurven[:,2]
        y_plot_3 = kurven[:,3]
        y_plot_4 = [self.h_gr - x*self.I_0 for x in x_plot]
        y_plot_5 = [self.h_N - x*self.I_0 for x in x_plot]

        plt.figure(dpi = 150)
        plt.plot(x_plot,y_plot_1,label = "Sohle")
        plt.plot(x_plot,y_plot_2,label = "Wasserspiegel")
        plt.plot(x_plot,y_plot_3,label = "Energielinie")
        plt.plot(x_plot,y_plot_4,'--',label="h_gr")
        plt.plot(x_plot,y_plot_5,'--',label="h_N")
        string_title = self.name + "-Spiegellinie: Sohle, Wasserspiegel, Energielinie"
        plt.title(string_title)
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.legend(prop={'size': 6})
        plt.show()

    def spiegelLinie(self):
        if (self.Fr_N<1):
            gefaelle_message = "Das Gefälle ist mild. Es können M1-Kurve, M2-Kurve oder M3-Kurve dargestellt werden."
            print(textwrap.fill(gefaelle_message,80))
            welche_Kurve = input("Wähle M1 mit der Taste '1' oder entsprechend M2/M3 mit '2' oder '3'. Wähle '4' für Beenden")
            welche_Kurve = int(welche_Kurve)
            while ((welche_Kurve == 1) or (welche_Kurve == 2) or (welche_Kurve == 3)):
                if(welche_Kurve == 1):
                    print("M1 benötigt einen Startwert größer als Normalabflusstiefe (hier:", round(self.h_N,3), "m):")
                    h_start = input("Bitte Startwert eingeben: ")
                    h_start = float(h_start)
                    self.__startWerte__(h_start,self.h_N)
                    self.__dxWerte__(-0.01,1.05)
                    self.name="M1"
                    self.berechnung()
                    welche_Kurve = input("Wähle M1 mit der Taste '1' oder entsprechend M2/M3 mit '2' oder '3'. Wähle '4' für Beenden")
                    welche_Kurve = int(welche_Kurve)
                elif(welche_Kurve == 2):
                    print("M2 benötigt einen Startwert kleiner als Normalabflusstiefe (hier:", round(self.h_N,3), "m) und größer als Grenzabflusstiefe (hier:", round(self.h_gr,3), "m):")
                    h_start = input("Bitte Startwert eingeben: ")
                    h_start = float(h_start)
                    self.__startWerte__(h_start,self.h_N)
                    self.__dxWerte__(-0.01,1.01)
                    self.name="M2"
                    self.berechnung()
                    welche_Kurve = input("Wähle M1 mit der Taste '1' oder entsprechend M2/M3 mit '2' oder '3'. Wähle '4' für Beenden")
                    welche_Kurve = int(welche_Kurve)
                elif(welche_Kurve == 3):
                    print("M3 benötigt einen Startwert kleiner als Grenzabflusstiefe (hier:", round(self.h_gr,3), "m):")
                    h_start = input("Bitte Startwert eingeben: ")
                    h_start = float(h_start)
                    self.__startWerte__(h_start,self.h_gr)
                    self.__dxWerte__(0.001,1.0)
                    self.name="M3"
                    self.berechnung()
                    welche_Kurve = input("Wähle M1 mit der Taste '1' oder entsprechend M2/M3 mit '2' oder '3'. Wähle '4' für Beenden")
                    welche_Kurve = int(welche_Kurve)
            print("ENDE")

        elif (self.Fr_N>1):
            gefaelle_message = "Das Gefälle ist steil. Es können S1-Kurve, S2-Kurve oder S3-Kurve dargestellt werden."
            print(textwrap.fill(gefaelle_message,80))
            welche_Kurve = input("Wähle S1 mit der Taste '1' oder entsprechend S2/S3 mit '2' oder '3'. Wähle '4' für Beenden")
            welche_Kurve = int(welche_Kurve)
            while ((welche_Kurve == 1) or (welche_Kurve == 2) or (welche_Kurve == 3)):
                if(welche_Kurve == 1):
                    print("S1 benötigt einen Startwert größer als Grenzabflusstiefe (hier:", round(self.h_gr,3), "m):")
                    h_start = input("Bitte Startwert eingeben: ")
                    h_start = float(h_start)
                    self.__startWerte__(h_start,self.h_gr)
                    self.__dxWerte__(-0.005,1.0001)
                    self.name="S1"
                    self.berechnung()
                    welche_Kurve = input("Wähle S1 mit der Taste '1' oder entsprechend S2/S3 mit '2' oder '3'. Wähle '4' für Beenden")
                    welche_Kurve = int(welche_Kurve)
                elif(welche_Kurve == 2):
                    print("S2 benötigt einen Startwert größer als Normalabflusstiefe (hier:", round(self.h_N,3), "m) und kleiner als Grenzabflusstiefe (hier:", round(self.h_gr,3), "m):")
                    h_start = input("Bitte Startwert eingeben: ")
                    h_start = float(h_start)
                    self.__startWerte__(h_start,self.h_N)
                    self.__dxWerte__(0.001,1.01)
                    self.name="S2"
                    self.berechnung()
                    welche_Kurve = input("Wähle S1 mit der Taste '1' oder entsprechend S2/S3 mit '2' oder '3'. Wähle '4' für Beenden")
                    welche_Kurve = int(welche_Kurve)
                elif(welche_Kurve == 3):
                    print("S3 benötigt einen Startwert kleiner als Normalabflusstiefe (hier:", round(self.h_N,3), "m):")
                    h_start = input("Bitte Startwert eingeben: ")
                    h_start = float(h_start)
                    self.__startWerte__(h_start,self.h_N)
                    self.__dxWerte__(0.05,1.001)
                    self.name="S3"
                    self.berechnung()
                    welche_Kurve = input("Wähle S1 mit der Taste '1' oder entsprechend S2/S3 mit '2' oder '3'. Wähle '4' für Beenden")
                    welche_Kurve = int(welche_Kurve)
            print("ENDE")

### Ende: Definition der Klasse Spiegellinie

if __name__ == "__main__":
    os.system('clear')

    start_message = "Dieses Programm berechnet für vorgegebenen Abfluss Q, Sohlgefälle I_0, Strickler-Wert k_str, Breite b und Böschungsneigung m die zugehörige Normalabflusstiefe, die Grenzabflusstiefe und die Froude-Zahl bei Normalabfluss. Außerdem wird der Gefälletyp (mild oder steil) bestimmt."
    print("\n")
    print(textwrap.fill(start_message,80))
    print("\n")
    start_message = "Darüberhinaus können Spiegellinien (M1, M2, M3, S1, S2, S3) berechnet und geplottet werden, je nach vorgegebenem Startwert einer Wassertiefe. Zu beachten ist dabei, dass bei strömendem Abfluss die zugrundeliegende Differentialgleichung von unterstrom nach oberstrom gelöst wird. Bei schießendem Abfluss erfolgt dies entsprechend andersherum. Dementsprechend werden die x-Achsen positiv oder negative Werte aufweisen bei den Plots."
    print(textwrap.fill(start_message,80))
    print("\n")

    start_message = "Beginne nun mit der Eingabe. Achte auf 'sinnvolle' Werte, sonst keine Garantie, dass das Programm durchläuft!"
    print(textwrap.fill(start_message,80))
    Q = input("Durchfluss Q in m^3/s: ")
    I_0 = input("Sohlgefälle I_0 [-]: ")
    k_str = input("Stricklerbeiwert k_str in m^(1/3) /s: ")
    b = input("Breite b in m: ")
    m = input("Böschungsneigung m [-]: ")
    Q=float(Q)
    b=float(b)
    m=float(m)
    k_str=float(k_str)
    I_0=float(I_0)
    sp = Spiegellinie(Q,b,I_0,k_str,m,"spiel")
    sp.normalAbfluss()
    sp.spiegelLinie()