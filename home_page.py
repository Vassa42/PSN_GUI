import PIL
from PIL import Image
import customtkinter
from contact_type import ContactType
from psn_evaluator import PSNEvaluator
from psn_preliminary import PSNPreliminary

customtkinter.set_appearance_mode("System")
customtkinter.set_default_color_theme("blue")

root = customtkinter.CTk()

# Carica l'immagine di background
image = PIL.Image.open("background.jpg")
background_image = customtkinter.CTkImage(image, size=(2500, 1500))

root.title('Home Page')
root.geometry('2500x1500')

def bg_resizer(e):
    if e.widget is root:
        i = customtkinter.CTkImage(image, size=(e.width, e.height))
        bg_label.configure(text="", image=i)

        # Ritaglia la parte dell'immagine di sfondo che corrisponde alla posizione della label
        update_label_bg()

def pass_to_contact_type():
    new_window = customtkinter.CTkToplevel(root)
    ContactType(new_window)

def pass_to_psn_preliminary():
    new_window = customtkinter.CTkToplevel(root)
    PSNPreliminary(new_window)

def pass_to_psn_evaluator():
    new_window = customtkinter.CTkToplevel(root)
    PSNEvaluator(new_window)

# Funzione per aggiornare lo sfondo della label in base alla sua posizione
def update_label_bg():
    # Posizione e dimensioni della label
    x = choose_module_label.winfo_x()
    y = choose_module_label.winfo_y()
    width = choose_module_label.winfo_width()
    height = choose_module_label.winfo_height()

    # Ritaglia la parte dell'immagine di background corrispondente alla label
    cropped_image = image.crop((x, y, x + width, y + height))

    # Converti l'immagine ritagliata in un oggetto utilizzabile da CTkImage
    cropped_ctk_image = customtkinter.CTkImage(cropped_image, size=(width, height))

    # Imposta l'immagine come sfondo della label
    choose_module_label.configure(image=cropped_ctk_image, text="Choose a Module:", text_color="red")

# Label per l'immagine di sfondo
bg_label = customtkinter.CTkLabel(root, text="", image=background_image)
bg_label.place(x=0, y=0)  # Posiziona lo sfondo nell'angolo in alto a sinistra

# Label con il testo "Choose a Module"
choose_module_label = customtkinter.CTkLabel(root, text="Choose a Module:", text_color="red", font=('Helvetica', 50))
choose_module_label.place(x=96, y=80)

# Pulsanti
contact_type_button = customtkinter.CTkButton(root, text="Contact Type", command=pass_to_contact_type, height=80, width=400, font=('Helvetica', 35), fg_color="transparent", hover_color="#a3a3a3")
contact_type_button.place(x=100, y=150)

psn_preliminary_button = customtkinter.CTkButton(root, text="PSN Preliminary", command=pass_to_psn_preliminary, height=80, width=400, font=('Helvetica', 35), fg_color="transparent", hover_color="#a3a3a3")
psn_preliminary_button.place(x=100, y=250)

psn_evaluator_button = customtkinter.CTkButton(root, text="PSN Evaluator", command=pass_to_psn_evaluator, height=80, width=400, font=('Helvetica', 35), fg_color="transparent", hover_color="#a3a3a3")
psn_evaluator_button.place(x=100, y=350)

# Bind per il ridimensionamento della finestra e aggiornamento dello sfondo della label
root.bind("<Configure>", bg_resizer)

# Aspetta che la finestra sia completamente inizializzata per ritagliare lo sfondo della label
root.after(100, update_label_bg)

root.mainloop()
