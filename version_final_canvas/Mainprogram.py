from pathlib import Path
from tkinter import Tk, Canvas, Entry, Text, Button, PhotoImage
import pandas as pd 
import datetime
import calendar
from PIL import Image, ImageTk
import numpy as np 
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import ttk
from tkinter import filedialog

def aide():
    w=Tk()
    w.title("Aide")
    w.rowconfigure(0, weight=1)
    w.columnconfigure(0, weight=1)
    w.configure(background="white")
    w.resizable(width=False, height=False)
    w.geometry("350x50+"+str(int(window.winfo_screenwidth()/2)-300)+"+"+str(int(window.winfo_screenheight()/2)-300))
    w.minsize(350,50)
    w.maxsize(350,50)
    Tiii=Label(w, text ="Verifiez que les colonnes de l'excel sont nomées :\n 'dateEmission' et 'heureEmission'.\n Verifiez qu'il sagit d'un fichier avec l'extension '.xlsx'.",bg='white')
    Tiii.place(x=0, y=0)
    w.mainloop()


def graph(event=None):
    
    if True :
        p=entry_1.get()
        def findDay(jour):
            born = datetime.datetime.strptime(jour, '%d %m %Y').weekday()
            return (calendar.day_name[born])
        def unique(list1):
            list_set = set(list1)
            unique_list = (list(list_set))
            return unique_list 

        df = pd.read_excel(p)
        D = df['dateEmission']
        H = df['heureEmission']
        heure = []
        date = []
        dateu =[]


        for i in H :
            k=str(i)
            heure.append(k[0:2])
        for i in D :          
            d = str(i.day)
            m = str(i.month)
            y = str(i.year)
            if len(str(i.day))== 1 :
                d = '0'+str(i.day)
            if len(str(i.month))== 1 :
                m = '0'+str(i.month)
            x = d+' '+m+' '+y
            date.append(findDay(x))
            
        
        UniDate = unique(D)

        for i in UniDate :          
            d = str(i.day)
            m = str(i.month)
            y = str(i.year)
            if len(str(i.day))== 1 :
                d = '0'+str(i.day)
            if len(str(i.month))== 1 :
                m = '0'+str(i.month)
            x = d+' '+m+' '+y
            dateu.append(findDay(x))
        lund =0
        mard =0
        mercred =0
        jeud =0
        vendred =0
        samed =0
        dimanch =0
        
        for i in dateu :
            if i == 'Monday':
                lund +=1
            elif i == 'Tuesday':
                mard +=1
            elif i == 'Wednesday':
                mercred +=1
            elif i == 'Thursday':
                jeud +=1
            elif i == 'Friday':
                vendred +=1
            elif i == 'Saturday':
                samed +=1
            elif i== 'Sunday':
                dimanch +=1
            
        lundi =[]
        mardi =[]
        mercredi =[]
        jeudi =[]
        vendredi =[]
        samedi =[]
        dimanche =[]

        mlundi =[]
        mmardi =[]
        mmercredi =[]
        mjeudi =[]
        mvendredi =[]
        msamedi =[]
        mdimanche =[]
        
        mdate=[]

        for i in range(len(date)):
            if date[i] == 'Monday':
                lundi.append(int(heure[i]))
            elif date[i] == 'Tuesday':
                mardi.append(int(heure[i]))
            elif date[i] == 'Wednesday':
                mercredi.append(int(heure[i]))
            elif date[i] == 'Thursday':
                jeudi.append(int(heure[i]))
            elif date[i] == 'Friday':
                vendredi.append(int(heure[i]))
            elif date[i] == 'Saturday':
                samedi.append(int(heure[i]))
            elif date[i] == 'Sunday':
                dimanche.append(int(heure[i]))
        for i in range (10,21):
            mlundi.append(lundi.count(i)/lund)
            mmardi.append(mardi.count(i)/mard)
            mmercredi.append(mercredi.count(i)/mercred)
            mjeudi.append(jeudi.count(i)/jeud)
            mvendredi.append(vendredi.count(i)/vendred)
            msamedi.append(samedi.count(i)/samed)
            mdimanche.append(dimanche.count(i)/dimanch)
        
        mdate.append(date.count('Monday')/lund)
        mdate.append(date.count('Tuesday')/mard)
        mdate.append(date.count('Wednesday')/mercred)
        mdate.append(date.count('Thursday')/jeud)
        mdate.append(date.count('Friday')/vendred)
        mdate.append(date.count('Saturday')/samed)
        mdate.append(date.count('Sunday')/dimanch)

            
    plt.close()
    q=entry_2.get()
    if (q=='Lundi'):
        x = np.arange(11)
        plt.bar(x, mlundi)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Lundi",loc='center',fontsize=8)
        plt.xlabel("Heures", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('10','11','12','13','14','15','16','17','18','19','20'))
        for x,y in zip(x,mlundi):
            plt.text(x, y-2, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()
    elif (q=='Mardi'):
        x = np.arange(11)
        plt.bar(x, mmardi)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Mardi",loc='center',fontsize=8)
        plt.xlabel("Heures", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('10','11','12','13','14','15','16','17','18','19','20'))
        for x,y in zip(x,mmardi):
            plt.text(x, y-2, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()
    elif (q=='Mercredi'):
        x = np.arange(11)
        plt.bar(x, mmercredi)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Mercredi",loc='center',fontsize=8)
        plt.xlabel("Heures", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('10','11','12','13','14','15','16','17','18','19','20'))
        for x,y in zip(x,mmercredi):
            plt.text(x, y-2, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()
    elif (q=='Jeudi'):
        x = np.arange(11)
        plt.bar(x, mjeudi)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Jeudi",loc='center',fontsize=8)
        plt.xlabel("Heures", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('10','11','12','13','14','15','16','17','18','19','20'))
        for x,y in zip(x,mjeudi):
            plt.text(x, y-2, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()
    elif (q=='Vendredi'):
        x = np.arange(11)
        plt.bar(x, mvendredi)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Vendredi",loc='center',fontsize=8)
        plt.xlabel("Heures", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('10','11','12','13','14','15','16','17','18','19','20'))
        for x,y in zip(x,mvendredi):
            plt.text(x, y-2, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()
    elif (q=='Samedi'):
        x = np.arange(11)
        plt.bar(x, msamedi)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Samedi",loc='center',fontsize=8)
        plt.xlabel("Heures", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('10','11','12','13','14','15','16','17','18','19','20'))
        for x,y in zip(x,msamedi):
            plt.text(x, y-2, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()
    elif (q=='Dimanche'):
        x = np.arange(11)
        plt.bar(x, mdimanche)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Dimanche",loc='center',fontsize=8)
        plt.xlabel("Heures", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('10','11','12','13','14','15','16','17','18','19','20'))
        for x,y in zip(x,mdimanche):
            plt.text(x, y-2, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()
    elif (q=='Semaine'):
        x = np.arange(7)
        plt.bar(x, mdate)
        plt.suptitle('Affluence des clients de la Détaxe', fontweight='bold', fontsize=18)
        plt.title("Semaine",loc='center',fontsize=8)
        plt.xlabel("Jour", fontsize=12)  
        plt.ylabel("Clients", fontsize=12)
        plt.xticks(x,('Lundi','Mardi','Mercredi','Jeudi','Vendredi','Samedi','Dimanche'))
        for x,y in zip(x,mdate):
            plt.text(x, y-8, '%.1f' % y, ha='center', va='bottom', fontsize=8, color='white', fontweight='bold')
        plt.show()

window = Tk()
v = StringVar()
fichier=['Lundi', 'Mardi' ,'Mercredi','Jeudi','Vendredi','Samedi','Dimanche','Semaine']
window.geometry("650x650")
window.configure(bg = "#FFFFFF")
window.title("Détaxe")

canvas = Canvas(window,bg = "#FFFFFF",height = 650,width = 650,bd = 0,highlightthickness = 0,relief = "ridge")

canvas.place(x = 0, y = 0)
canvas.create_text(11.0,16.0,anchor="nw",text="Affluence App",fill="#000000",font=("Inter", 16 * -1))
def search():
    file = filedialog.askopenfile(parent=window,mode='rb',title="Selectionnez l'excel")
    d=str(file)
    p=d[26:]
    a=p[:-2]
    v.set(a)
button_image_1 = PhotoImage(file="button_1.png")
button_1 = Button(image=button_image_1,borderwidth=0,highlightthickness=0,command=search,relief="flat")
button_1.place(x=482.0,y=462.0,width=143.0,height=28.0)

button_image_2 = PhotoImage(file="button_2.png")
button_2 = Button(image=button_image_2,borderwidth=0,highlightthickness=0,command=graph,relief="flat")
button_2.place(x=482.0,y=536.0,width=143.0,height=28.0)

entry_image_1 = PhotoImage(file="entry_1.png")
entry_bg_1 = canvas.create_image(325.0,475.5,image=entry_image_1)
entry_1 = Entry(bd=0,bg="#F1F1F1",highlightthickness=0,textvariable = v )
entry_1.place(x=210.0,y=476.0,width=230.0,height=20.0)

canvas.create_text(204.0,452.0,anchor="nw",text="Nom du fichier",fill="#000000",font=("Inter Light", 14 * -1))

entry_image_2 = PhotoImage(file="entry_2.png")
entry_bg_2 = canvas.create_image(325.0,549.5,image=entry_image_2)
entry_2 = ttk.Combobox(values=fichier)
entry_2.place(x=210.0,y=550.0,width=230.0,height=20.0)

canvas.create_text(204.0,524.0,anchor="nw",text="Jour",fill="#000000",font=("Inter Light", 14 * -1))

image_image_1 = PhotoImage(file="image_1.png")
image_1 = canvas.create_image(325.0,249.0,image=image_image_1)

canvas.create_rectangle(11.0,36.0,374.0,36.0,fill="#000000",outline="")

button_image_3 = PhotoImage(file="button_3.png")
button_3 = Button(image=button_image_3,bg="white",borderwidth=0,highlightthickness=0,command=aide,relief="flat")
button_3.place(x=61.0,y=529.0,width=33.0,height=44.0)
window.resizable(False, False)
window.mainloop()
