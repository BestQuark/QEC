from tkinter import *
import syndrome as syn
import os

"""
GUI that uses syndrome.py to calculate errors associated to error syndromes of different stabilizer codes
Author: @BestQuark
"""

os.system('clear')

def calculateErrors():
    txt = syn.syndromes_symplectic2txt(syn.minimum_weight_error(syn.classify_errors(syndromeTextbox.get())))
    
    titleSyn = Label(root, text='Error Syndrome')
    titleErr = Label(root, text='Pauli errors')
    titleSyn.grid(row=6, column=0)
    titleErr.grid(row=6, column=1)
    titleSyn.configure(font=("Times New Roman", 12, "bold"), bg='white')
    titleErr.configure(font=("Times New Roman", 12, "bold"), bg='white')
    
    syndrosList = Listbox(root, width = 20)
    syndrosList.grid(row=7, column=0)
    syndrosList.configure(font=("Times New Roman", 12))
    
    errsList = Listbox(root, width = 50)
    errsList.grid(row=7, column=1)
    errsList.configure(font=("Times New Roman", 12))
    
    for indx, err in enumerate(txt.keys()):
        syndrosList.insert(END, str(indx+1)+')     '+str(err))
        
        if type(txt[err]) is str:
            ccc = str(txt[err])
        else:
            ccc = " - ".join(txt[err])
        errsList.insert(END, str(indx+1)+')     '+ ccc)
        
def exampleSelected(event):
    if exampleClicked.get()==stabExamples[0]:
        syndromeTextbox.delete(0, END)
        syndromeTextbox.insert(0, "XXXXIII;XXIXXII;XIXIXII;ZZZZIII; ZZIZZII; ZIZIZII")
        
    elif exampleClicked.get()==stabExamples[1]:
        syndromeTextbox.delete(0, END)
        syndromeTextbox.insert(0, "IZXXZ; ZIZXX; XZIZX; XXZIZ")
        
    elif exampleClicked.get()==stabExamples[2]:
        syndromeTextbox.delete(0, END)
        syndromeTextbox.insert(0, "XXI;IXX")  
    
    exampleClicked.set("Examples")

    
root = Tk()

root.title("Syndrome QEC")
root.iconbitmap("iconAtom.ico")
root.geometry('700x450')
root.configure(background='white')


stabExamples = ["Steane code", "5 qubit code", "Repetition code"]

exampleClicked = StringVar()
exampleClicked.set("Examples")

drop = OptionMenu(root, exampleClicked, *stabExamples, command=exampleSelected)
drop.configure(font=("Times New Roman", 12), bg='blue', fg='white')
drop.grid(row=5, column=2, sticky="ew")

appTitle = Label(root, text='Syndrome Calculator')
appTitle.grid(row=0, column=1, pady=10)
appTitle.configure(font=("Times New Roman", 25, "bold"), bg='white')

instructionsLb = Label(root, text = 'To calculate Pauli errors associated with an error syndrome, enter the stabilizer of the code as \n "S1; S2; S3; ... ; Sn" where Si is a product of Pauli matrices. ')
instructionsLb.grid(row=1, column=0, columnspan = 2, rowspan=3, pady=15)
instructionsLb.configure(font=("Times New Roman", 12), bg='white')

syndromeLabel = Label(root, text='Stabilizer:')
syndromeLabel.grid(row=4, column=0, pady=10)
syndromeLabel.configure(font=("Times New Roman", 15, "bold") , bg='white')

#syndromeLabel.pack()

syndromeTextbox = Entry(root, width=50,  borderwidth=2)
syndromeTextbox.grid(row=4, column=1)
syndromeTextbox.configure(font=("Times New Roman", 12))
#syndromeTextbox.pack()

errorButton = Button(root, text='Calculate', command=calculateErrors, bg='blue', fg='white')
errorButton.grid(row=4, column=2)
errorButton.configure(font=("Times New Roman", 12))
#errorButton.pack()

root.mainloop()