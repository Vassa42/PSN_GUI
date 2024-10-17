import PIL
from PIL import Image
import customtkinter
from contact_type import ContactType
from psn_evaluator import PSNEvaluator
from psn_preliminary import PSNPreliminary

# Configurazione iniziale dell'interfaccia grafica
customtkinter.set_appearance_mode("System")  # Modalità chiara/scura automatica
customtkinter.set_default_color_theme("blue")

root = customtkinter.CTk()

# Carica l'immagine di background
image = PIL.Image.open("background.jpg")
background_image = customtkinter.CTkImage(image, size=(2500, 1500))

root.title('Home Page')
root.geometry('2500x1500')

# Funzione per ridimensionare l'immagine di sfondo quando la finestra cambia dimensione
def bg_resizer(e):
    if e.widget is root:
        i = customtkinter.CTkImage(image, size=(e.width, e.height))  # Aggiorna immagine di sfondo
        bg_label.configure(text="", image=i)
        update_label_bg()  # Aggiorna la parte ritagliata dello sfondo per la label

# Funzioni per aprire nuove finestre relative ai moduli
def pass_to_contact_type():
    new_window = customtkinter.CTkToplevel(root)
    ContactType(new_window)

def pass_to_psn_preliminary():
    new_window = customtkinter.CTkToplevel(root)
    PSNPreliminary(new_window)

def pass_to_psn_evaluator():
    new_window = customtkinter.CTkToplevel(root)
    PSNEvaluator(new_window)

# Funzione per aggiornare lo sfondo della label in base alla sua posizione nella finestra
def update_label_bg():
    # Posizione e dimensioni della label "Choose a Module"
    x = choose_module_label.winfo_x()
    y = choose_module_label.winfo_y()
    width = choose_module_label.winfo_width()
    height = choose_module_label.winfo_height()

    # Ritaglia la parte dell'immagine di sfondo corrispondente alla label
    cropped_image = image.crop((x, y, x + width, y + height))
    cropped_ctk_image = customtkinter.CTkImage(cropped_image, size=(width, height))

    # Imposta l'immagine come sfondo della label
    choose_module_label.configure(image=cropped_ctk_image, text="Choose a Module:", text_color="red")

# Imposta il background della finestra principale
bg_label = customtkinter.CTkLabel(root, text="", image=background_image)
bg_label.place(x=0, y=0)

# Label con il testo "Choose a Module"
choose_module_label = customtkinter.CTkLabel(root, text="Choose a Module:", text_color="red", font=('Helvetica', 50))
choose_module_label.place(x=96, y=80)

# Pulsanti per accedere ai diversi moduli
contact_type_button = customtkinter.CTkButton(root, text="Contact Type", command=pass_to_contact_type, height=80, width=400, font=('Helvetica', 35), fg_color="transparent", hover_color="#a3a3a3")
contact_type_button.place(x=100, y=150)

psn_preliminary_button = customtkinter.CTkButton(root, text="PSN Preliminary", command=pass_to_psn_preliminary, height=80, width=400, font=('Helvetica', 35), fg_color="transparent", hover_color="#a3a3a3")
psn_preliminary_button.place(x=100, y=250)

psn_evaluator_button = customtkinter.CTkButton(root, text="PSN Evaluator", command=pass_to_psn_evaluator, height=80, width=400, font=('Helvetica', 35), fg_color="transparent", hover_color="#a3a3a3")
psn_evaluator_button.place(x=100, y=350)

# Aggiorna l'immagine di sfondo e la label quando la finestra viene ridimensionata
root.bind("<Configure>", bg_resizer)

# Esegui l'aggiornamento dello sfondo della label dopo che la finestra è stata inizializzata
root.after(100, update_label_bg)

root.mainloop()
