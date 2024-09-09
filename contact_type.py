import customtkinter
from customtkinter import filedialog
from contact_calculator_sb import ContactCalculatorSb
from contact_calculator_hb import ContactCalculatorHb
from contact_calculator_hc import ContactCalculatorHc


class ContactType:

    def __init__(self, root):
        self.root = root
        self.root.title("Contact type")
        self.root.configure(bg='LightSkyBlue1')

        self.outer_frame = customtkinter.CTkFrame(self.root, fg_color='LightSkyBlue1', corner_radius=10)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        self.frame = customtkinter.CTkFrame(self.outer_frame, fg_color='LightSkyBlue1', corner_radius=10)
        self.frame.pack(expand=True, padx=10, pady=10)

        self.label_output = customtkinter.CTkLabel(self.frame, text="Select the contact's type:", font=("Helvetica", 45), text_color="#000000")
        self.label_output.grid(row=0, column=0, padx=10, pady=10)

        # Variabili per i checkbox
        self.var_salt_bridges = customtkinter.IntVar()
        self.var_hydrophobic_clusters = customtkinter.IntVar()
        self.var_hydrogen_bonds = customtkinter.IntVar()

        # Checkbox per i vari moduli
        self.check_SB = customtkinter.CTkCheckBox(self.frame, text="Salt Bridges", variable=self.var_salt_bridges,
                                                  font=("Helvetica", 35), text_color="#000000", command=self.update_button_state)
        self.check_HC = customtkinter.CTkCheckBox(self.frame, text="Hydrophobic Clusters", variable=self.var_hydrophobic_clusters,
                                                  font=("Helvetica", 35), text_color="#000000", command=self.update_button_state)
        self.check_HB = customtkinter.CTkCheckBox(self.frame, text="Hydrogen Bonds", variable=self.var_hydrogen_bonds,
                                                  font=("Helvetica", 35), text_color="#000000", command=self.update_button_state)

        self.check_SB.grid(row=1, column=0, sticky='w', padx=10)
        self.check_HC.grid(row=2, column=0, sticky='w', padx=10)
        self.check_HB.grid(row=3, column=0, sticky='w', padx=10)

        # Bottone per selezionare i file e avviare i calcolatori
        self.button_select = customtkinter.CTkButton(self.frame, text="Select", command=self.open_file_selectors,
                                                     state="disabled", font=("Helvetica", 35),
                                                     fg_color="#333333",  # Dark color for button
                                                     hover_color="#555555")  # Slightly lighter hover color
        self.button_select.grid(row=4, column=0, columnspan=3, pady=20)

        self.label_result = customtkinter.CTkLabel(self.frame, text="", font=("Helvetica", 35), fg_color='LightSkyBlue1', text_color="#000000")
        self.label_result.grid(row=5, column=0, columnspan=3, pady=10)

        # Aggiungi campi per memorizzare file
        self.topology_file = None
        self.trajectory_file = None
        self.output_directory = None

    # Aggiorna lo stato del bottone Select
    def update_button_state(self):
        if self.var_salt_bridges.get() or self.var_hydrophobic_clusters.get() or self.var_hydrogen_bonds.get():
            self.button_select.configure(state="normal")
        else:
            self.button_select.configure(state="disabled")

    # Funzione per aprire i file di topologia, traiettoria e la cartella di output
    def open_file_selectors(self):
        self.topology_file = filedialog.askopenfilename(title="Select the topology file",
                                                        filetypes=[("PDB files", "*.pdb"), ("GRO files", "*.gro")])
        self.trajectory_file = filedialog.askopenfilename(title="Select the trajectory file",
                                                          filetypes=[("PDB files", "*.pdb"), ("XTC files", "*.xtc"),
                                                                     ("DCD files", "*.dcd")])
        self.output_directory = filedialog.askdirectory(title="Select Output Directory")

        # Apri i calcolatori selezionati
        if self.topology_file and self.trajectory_file and self.output_directory:
            self.open_contact_calculators()

    # Funzione per aprire le finestre dei calcolatori selezionati
    def open_contact_calculators(self):
        if self.var_salt_bridges.get() and self.var_hydrophobic_clusters.get() and self.var_hydrogen_bonds.get():
            self.open_calculator(ContactCalculatorSb)
            self.open_calculator(ContactCalculatorHc)
            self.open_calculator(ContactCalculatorHb)

        elif self.var_salt_bridges.get() and self.var_hydrophobic_clusters.get():
            self.open_calculator(ContactCalculatorSb)
            self.open_calculator(ContactCalculatorHc)

        elif self.var_salt_bridges.get() and self.var_hydrogen_bonds.get():
            self.open_calculator(ContactCalculatorSb)
            self.open_calculator(ContactCalculatorHb)

        elif self.var_hydrophobic_clusters.get() and self.var_hydrogen_bonds.get():
            self.open_calculator(ContactCalculatorHc)
            self.open_calculator(ContactCalculatorHb)

        elif self.var_salt_bridges.get():
            self.open_calculator(ContactCalculatorSb)

        elif self.var_hydrophobic_clusters.get():
            self.open_calculator(ContactCalculatorHc)

        elif self.var_hydrogen_bonds.get():
            self.open_calculator(ContactCalculatorHb)

    # Funzione per aprire una nuova finestra del calcolatore
    def open_calculator(self, calculator_class):
        new_window = customtkinter.CTkToplevel(self.root)
        calculator = calculator_class(new_window, self.topology_file, self.trajectory_file, self.output_directory)


if __name__ == "__main__":
    customtkinter.set_appearance_mode("dark")  # Modalità scura
    customtkinter.set_default_color_theme("green")  # Tema verde

    root = customtkinter.CTk()
    app = ContactType(root)
    root.mainloop()
