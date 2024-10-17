import customtkinter
import customtkinter as ctk

class ContactCalculatorHc:
    def __init__(self, root, topology_file, trajectory_file, output_directory):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.output_directory = output_directory
        self.root = root
        self.root.title("Contact calculator Hydrophobic Clusters")

        # Frame esterno che contiene tutti gli elementi dell'interfaccia, con un colore personalizzato
        self.outer_frame = ctk.CTkFrame(root, fg_color="DeepSkyBlue2", corner_radius=10)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per posizionare i widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="DeepSkyBlue2")
        self.frame.pack(expand=True, padx=10, pady=10)

        # Etichetta e campo di input per inserire la distanza di cut-off
        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("Helvetica", 40), text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10)

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("Helvetica", 30))
        self.entry_cutOff.insert(0, "5.0")  # Valore predefinito per il cut-off
        self.entry_cutOff.grid(row=1, column=1, padx=10, pady=10)

        # Bottone per confermare il valore del cut-off
        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30), command=self.confirm_cutoff, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=1, column=2, padx=10, pady=10)

        # Etichetta per i residui idrofobici
        self.label_residues = ctk.CTkLabel(self.frame, text="Hydrophobic residues:", font=("Helvetica", 40), text_color="black")
        self.label_residues.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        # Variabili per tenere traccia dei checkbox selezionati per i residui idrofobici
        self.var_ala = customtkinter.IntVar()
        self.var_val = customtkinter.IntVar()
        self.var_leu = customtkinter.IntVar()
        self.var_ile = customtkinter.IntVar()
        self.var_phe = customtkinter.IntVar()
        self.var_pro = customtkinter.IntVar()
        self.var_trp = customtkinter.IntVar()
        self.var_met = customtkinter.IntVar()
        self.var_tyr = customtkinter.IntVar()
        self.var_gly = customtkinter.IntVar()

        # Checkbox per la selezione dei residui idrofobici
        self.check_ala = customtkinter.CTkCheckBox(self.frame, text="ALA", variable=self.var_ala, font=("Helvetica", 30), text_color="#000000")
        self.check_val = customtkinter.CTkCheckBox(self.frame, text="VAL", variable=self.var_val, font=("Helvetica", 30), text_color="#000000")
        self.check_leu = customtkinter.CTkCheckBox(self.frame, text="LEU", variable=self.var_leu, font=("Helvetica", 30), text_color="#000000")
        self.check_ile = customtkinter.CTkCheckBox(self.frame, text="ILE", variable=self.var_ile, font=("Helvetica", 30), text_color="#000000")
        self.check_phe = customtkinter.CTkCheckBox(self.frame, text="PHE", variable=self.var_phe, font=("Helvetica", 30), text_color="#000000")
        self.check_pro = customtkinter.CTkCheckBox(self.frame, text="PRO", variable=self.var_pro, font=("Helvetica", 30), text_color="#000000")
        self.check_trp = customtkinter.CTkCheckBox(self.frame, text="TRP", variable=self.var_trp, font=("Helvetica", 30), text_color="#000000")
        self.check_met = customtkinter.CTkCheckBox(self.frame, text="MET", variable=self.var_met, font=("Helvetica", 30), text_color="#000000")
        self.check_tyr = customtkinter.CTkCheckBox(self.frame, text="TYR", variable=self.var_tyr, font=("Helvetica", 30), text_color="#000000")
        self.check_gly = customtkinter.CTkCheckBox(self.frame, text="GLY", variable=self.var_gly, font=("Helvetica", 30), text_color="#000000")

        # Posizionamento dei checkbox per i residui
        self.check_ala.grid(row=3, column=0, sticky='w', padx=10)
        self.check_val.grid(row=4, column=0, sticky='w', padx=10)
        self.check_leu.grid(row=5, column=0, sticky='w', padx=10)
        self.check_ile.grid(row=6, column=0, sticky='w', padx=10)
        self.check_phe.grid(row=7, column=0, sticky='w', padx=10)
        self.check_pro.grid(row=8, column=0, sticky='w', padx=10)
        self.check_trp.grid(row=9, column=0, sticky='w', padx=10)
        self.check_met.grid(row=10, column=0, sticky='w', padx=10)
        self.check_tyr.grid(row=11, column=0, sticky='w', padx=10)
        self.check_gly.grid(row=12, column=0, sticky='w', padx=10)

        # Etichetta per la selezione del file delle masse FF (forze di campo)
        self.label_ff_masses = ctk.CTkLabel(self.frame, text="FF_masses file:", font=("Helvetica", 40), text_color="black")
        self.label_ff_masses.grid(row=2, column=1, padx=10, pady=10, sticky='w')

        # Lista di masse di forze di campo disponibili
        self.ff_masses = ["charmm27", "amber03", "amber94", "amber96", "amber99", "amber99sb", "amber99as-iidn", "amberGS", "encads", "encadv", "gromos43a1", "gromos43a2", "gromos45a3", "gromos53a5", "groms53a6", "oplsaa"]
        self.selected_ff_mass = ctk.StringVar()

        # Menu a tendina per la selezione delle masse FF con font personalizzato
        self.combobox_class = ctk.CTkComboBox(self.frame, values=self.ff_masses, variable=self.selected_ff_mass, state='readonly', width=200, font=("helvetica", 30), dropdown_font=("helvetica", 30))
        self.combobox_class.grid(row=2, column=2, padx=10, pady=10)
        self.combobox_class.set(self.ff_masses[0])  # Imposta valore di default

        # Bottone per avviare il processo di calcolo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), fg_color="gray30", width=200)
        self.button_process.grid(row=13, column=0, columnspan=3, pady=20)

        # Etichetta per mostrare il risultato dell'elaborazione
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40), text_color="black")
        self.label_result.grid(row=14, column=0, columnspan=3, pady=10)

    # Funzione per ottenere il valore del cut-off inserito dall'utente
    def confirm_cutoff(self):
        self.entry_cutOff.get()


# Avvio della GUI
if __name__ == "__main__":
    ctk.set_appearance_mode("dark")  # Modalità scura
    ctk.set_default_color_theme("blue")  # Imposta un tema di default

    root = ctk.CTk()  # Crea la finestra principale con customtkinter
    app = ContactCalculatorHc(root, "topology.pdb", "trajectory.dcd", "output_directory")
    root.mainloop()
