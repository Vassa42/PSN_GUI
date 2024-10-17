import customtkinter as ctk

class ContactCalculatorHb:
    def __init__(self, root, topology_file, trajectory_file, output_directory):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.output_directory = output_directory
        self.root = root
        self.root.title("Contact Calculator Hydrogen Bonds")  # Titolo della finestra principale

        # Imposta tema scuro di default
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")  # Imposta un tema colore blu

        # Frame esterno che contiene l'interfaccia con bordi e colore di sfondo personalizzato
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="DodgerBlue2")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per allineare i vari widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="DodgerBlue2")
        self.frame.pack(expand=True, padx=10, pady=10)

        # Label e campo di input per la distanza di cut-off
        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("helvetica", 40), text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("helvetica", 30))
        self.entry_cutOff.insert(0, "3.5")  # Valore di default per la distanza di cut-off
        self.entry_cutOff.grid(row=1, column=1, padx=10, pady=10)

        # Bottone per confermare il valore di cut-off
        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30), command=self.confirm_cutoff, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=1, column=2, padx=10, pady=10)

        # Label e campo di input per l'angolo del legame idrogeno
        self.label_angle = ctk.CTkLabel(self.frame, text="Insert the bond's angle (D):", font=("helvetica", 40), text_color="black")
        self.label_angle.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        self.entry_angle = ctk.CTkEntry(self.frame, width=200, justify="center", font=("helvetica", 30))
        self.entry_angle.insert(0, "120.0")  # Valore di default per l'angolo del legame
        self.entry_angle.grid(row=2, column=1, padx=10, pady=10)

        # Bottone per confermare il valore dell'angolo
        self.button_angle = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30), command=self.confirm_angle, fg_color="gray30", width=200)
        self.button_angle.grid(row=2, column=2, padx=10, pady=10)

        # Label per selezionare la classe di legami idrogeno da analizzare
        self.label_class = ctk.CTkLabel(self.frame, text="Select the class to analyze:", font=("Helvetica", 40), text_color="black")
        self.label_class.grid(row=3, column=0, padx=10, pady=10, sticky='w')

        # Definisce le opzioni per la classe di legami idrogeno da analizzare
        self.classes = ["All", "Mc-Mc", "Mc-Sc", "Sc-Sc"]
        self.selected_class = ctk.StringVar()

        # Menu a tendina per selezionare la classe di legami, con font personalizzato
        self.combobox_class = ctk.CTkComboBox(self.frame, values=self.classes, variable=self.selected_class, state='readonly', width=200, font=("helvetica", 30), dropdown_font=("helvetica", 30))
        self.combobox_class.grid(row=3, column=1, padx=10, pady=10)
        self.combobox_class.set(self.classes[0])  # Imposta il valore di default

        # Label e campo per la selezione del file di legami idrogeno
        self.label_hb = ctk.CTkLabel(self.frame, text="Select the hydrogen bonds' file:", font=("Helvetica", 40), text_color="black")
        self.label_hb.grid(row=4, column=0, padx=10, pady=10)

        # Campo di input per inserire il percorso del file dei legami idrogeno
        self.entry_hb = ctk.CTkEntry(self.frame, width=500, justify="center", font=("Helvetica", 25))
        self.entry_hb.insert(0, "/home/vassa42/Documenti/Tirocinio-Tesi/PSN_GUI/hydrogen_bonds.ini")  # Valore di default del file
        self.entry_hb.grid(row=4, column=1, padx=10, pady=10)

        # Bottone per navigare e selezionare un file INI
        self.button_hb = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 30), command=self.browse_hb_file, width=200, height=40, fg_color="gray30")
        self.button_hb.grid(row=4, column=2, padx=10, pady=10)

        # Bottone per avviare il processo di calcolo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), fg_color="gray30", width=200)
        self.button_process.grid(row=5, column=0, columnspan=3, pady=20)

        # Label per visualizzare i risultati o eventuali messaggi di errore
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40))
        self.label_result.grid(row=6, column=0, columnspan=3, pady=10)

    # Metodo per ottenere il valore della distanza di cut-off
    def confirm_cutoff(self):
        self.entry_cutOff.get()

    # Metodo per ottenere il valore dell'angolo
    def confirm_angle(self):
        self.entry_angle.get()

    # Metodo per navigare e selezionare il file dei legami idrogeno
    def browse_hb_file(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the hydrogen bonds' file", filetypes=[("INI files", "*.ini"), ("All files", "*.*")])
        if filename:
            self.entry_hb.delete(0, ctk.END)
            self.entry_hb.insert(0, filename)  # Aggiorna il campo di input con il percorso del file selezionato

# Avvio dell'applicazione GUI
if __name__ == "__main__":
    root = ctk.CTk()
    app = ContactCalculatorHb(root, "topology.psf", "trajectory.dcd", "output_directory")
    root.mainloop()
