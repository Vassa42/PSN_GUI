import customtkinter as ctk

class ContactCalculatorSb:
    def __init__(self, root, topology_file, trajectory_file, output_directory):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.output_directory = output_directory
        self.root = root
        self.root.title("Contact Calculator: Salt Bridges")  # Imposta il titolo della finestra

        # Crea un frame esterno per contenere gli elementi dell'interfaccia, con colori personalizzati
        self.outer_frame = ctk.CTkFrame(root, fg_color='SteelBlue2', corner_radius=10)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per aggiungere i widget dell'interfaccia
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color='SteelBlue2')
        self.frame.pack(expand=True, padx=10, pady=10)

        # Etichetta e campo di input per inserire la distanza di cut-off (in Ångström)
        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("Helvetica", 40), fg_color='SteelBlue2', text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("Helvetica", 30))
        self.entry_cutOff.insert(0, "4.5")  # Valore di default per la distanza di cut-off
        self.entry_cutOff.grid(row=1, column=1, padx=10, pady=10)

        # Bottone per confermare la distanza di cut-off inserita
        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30), command=self.confirm_cutoff, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=1, column=2, padx=10, pady=10)

        # Etichetta e campo di input per selezionare il file dei gruppi caricati
        self.label_charged = ctk.CTkLabel(self.frame, text="Select the charged groups' file:", font=("Helvetica", 40), text_color="black")
        self.label_charged.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        self.entry_charged = ctk.CTkEntry(self.frame, width=500, justify="center", font=("Helvetica", 25))
        # Inserisce un percorso di file predefinito come placeholder
        self.entry_charged.insert(0, "/home/vassa42/Documenti/Tirocinio-Tesi/PSN_GUI/charged_groups.ini")
        self.entry_charged.grid(row=2, column=1, padx=10, pady=10)

        # Bottone per sfogliare e selezionare il file dei gruppi caricati
        self.button_charged = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 30), command=self.browse_charged_groups_file, width=200, height=40, fg_color="gray30")
        self.button_charged.grid(row=2, column=2, padx=10, pady=10)

        # Bottone per avviare il processo di calcolo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), fg_color="gray30", width=200)
        self.button_process.grid(row=3, column=0, columnspan=3, pady=20)

        # Etichetta per visualizzare i risultati o eventuali errori durante l'elaborazione
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40), fg_color='SteelBlue2', text_color="black")
        self.label_result.grid(row=4, column=0, columnspan=3, pady=10)

    # Funzione per confermare e ottenere il valore del cut-off inserito
    def confirm_cutoff(self):
        self.entry_cutOff.get()

    # Funzione per aprire una finestra di dialogo e selezionare il file dei gruppi caricati
    def browse_charged_groups_file(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the charged groups' file", filetypes=[("INI files", "*.ini"), ("All files", "*.*")])
        if filename:
            self.entry_charged.delete(0, ctk.END)  # Cancella il valore attuale del campo
            self.entry_charged.insert(0, filename)  # Inserisce il percorso del file selezionato


# Avvio della GUI principale
if __name__ == "__main__":
    ctk.set_appearance_mode("dark")  # Imposta la modalità scura
    ctk.set_default_color_theme("blue")  # Imposta il tema blu di default

    root = ctk.CTk()
    app = ContactCalculatorSb(root, "topology.pdb", "trajectory.dcd", "output_directory")
    root.mainloop()
