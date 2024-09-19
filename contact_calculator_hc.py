import customtkinter
import customtkinter as ctk
import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree
import os

class ContactCalculatorHc:
    def __init__(self, root, topology_file, trajectory_file, output_directory):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.output_directory = output_directory
        self.root = root
        self.root.title("Contact calculator Hydrophobic Clusters")

        # Configurazione della finestra principale in modalità scura
        self.outer_frame = ctk.CTkFrame(root, fg_color="DeepSkyBlue2", corner_radius=10)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="DeepSkyBlue2")
        self.frame.pack(expand=True, padx=10, pady=10)

        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("Helvetica", 40), text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10)

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("Helvetica", 30))
        self.entry_cutOff.insert(0, "5.0")
        self.entry_cutOff.grid(row=1, column=1, padx=10, pady=10)

        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30), command=self.confirm_cutoff, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=1, column=2, padx=10, pady=10)

        self.label_residues = ctk.CTkLabel(self.frame, text="Hydrophobic residues:", font=("Helvetica", 40), text_color="black")
        self.label_residues.grid(row=2, column=0, padx=10, pady=10, sticky='w')

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

        self.label_ff_masses = ctk.CTkLabel(self.frame, text="FF_masses file:", font=("Helvetica", 40), text_color="black")
        self.label_ff_masses.grid(row=2, column=1, padx=10, pady=10, sticky='w')

        self.ff_masses = ["charmm27", "amber03", "amber94", "amber96", "amber99", "amber99sb", "amber99as-iidn", "amberGS", "encads", "encadv", "gromos43a1", "gromos43a2", "gromos45a3", "gromos53a5", "groms53a6", "oplsaa"]
        self.selected_ff_mass = ctk.StringVar()

        # Menu a tendina con il font personalizzato
        self.combobox_class = ctk.CTkComboBox(self.frame, values=self.ff_masses, variable=self.selected_ff_mass, state='readonly', width=200, font=("helvetica", 30), dropdown_font=("helvetica", 30))
        self.combobox_class.grid(row=2, column=2, padx=10, pady=10)
        self.combobox_class.set(self.ff_masses[0])

        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), fg_color="gray30", width=200)
        self.button_process.grid(row=13, column=0, columnspan=3, pady=20)

        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40), text_color="black")
        self.label_result.grid(row=14, column=0, columnspan=3, pady=10)

    def confirm_cutoff(self):
        self.entry_cutOff.get()
"""        
    def preprocess_file(self):
        topology_file = self.topology_file
        trajectory_file = self.trajectory_file
        cutoff_distance = float(self.entry_cutOff.get())
        output_directory = self.output_directory

        if topology_file and trajectory_file and cutoff_distance:
            self.process_file(topology_file, trajectory_file, cutoff_distance, output_directory)
            self.label_result.configure(text="Processing completed successfully!", text_color="white")
        else:
            self.label_result.configure(text="Please fill all fields correctly.", text_color="red")

    def process_file(self, topology_file, trajectory_file, cutoff_distance, output_directory):
        # Carica la topologia e la traiettoria usando MDAnalysis
        u = mda.Universe(topology_file, trajectory_file)

        # Ottieni il nome del file di topologia senza l'estensione
        file_prefix = os.path.splitext(os.path.basename(topology_file))[0]

        # Seleziona i residui idrofobici
        hydrophobic_residues = u.select_atoms("resname ALA VAL LEU ILE MET PHE TRP PRO")

        # Dizionario per memorizzare la persistenza delle interazioni
        interaction_dict = {}
        frame_count = 0

        # Trova le coppie di residui che interagiscono entro il cutoff per ogni frame
        for ts in u.trajectory:
            frame_count += 1
            # Crea KD-tree per una ricerca spaziale efficiente
            kd_tree = cKDTree(hydrophobic_residues.positions)

            # Trova coppie di residui entro il raggio di cutoff
            pairs = kd_tree.query_pairs(r=cutoff_distance)

            # Aggiungi o aggiorna le coppie che formano contatti idrofobici
            for i, j in pairs:
                # Usa una tupla ordinata per evitare duplicati (i,j) e (j,i)
                pair = tuple(sorted((i, j)))
                if pair not in interaction_dict:
                    interaction_dict[pair] = 1  # Nuova interazione osservata
                else:
                    interaction_dict[pair] += 1  # Incrementa la persistenza

        # Numero totale di interazioni effettivamente osservate
        n_interactions = len(interaction_dict)
        print(f"Numero di interazioni idrofobiche uniche osservate: {n_interactions}")

        # Crea una matrice di persistenza di dimensioni n_interactions x n_interactions
        persistence_matrix = np.zeros((len(hydrophobic_residues), len(hydrophobic_residues)))

        # Popola la matrice di persistenza con i dati raccolti
        for (i, j), persistence in interaction_dict.items():
            persistence_percentage = (persistence / frame_count) * 100
            persistence_matrix[i, j] = persistence_percentage
            persistence_matrix[j, i] = persistence_percentage  # La matrice è simmetrica

        # Salva i contatti idrofobici in un file di testo con lo stesso formato dei ponti salini
        contacts_file = os.path.join(output_directory, f"{file_prefix}_contacts_HC.txt")
        persistence_matrix_file = os.path.join(output_directory, f"{file_prefix}_matrix_HC.dat")

        with open(contacts_file, 'w') as f:
            f.write("#CHAIN\tRES1_ID\tRES1_NAME\tRES1_GROUP\tRES2_ID\tRES2_NAME\tRES2_GROUP\tOCCURRENCE_PERC\n")
            for (i, j), persistence in interaction_dict.items():
                # Solo contatti con persistenza
                res1 = hydrophobic_residues[i]
                res2 = hydrophobic_residues[j]

                chain1 = res1.segid
                chain2 = res2.segid
                res1_id = res1.resid
                res1_name = res1.resname
                res2_id = res2.resid
                res2_name = res2.resname

                # Descrizione del gruppo laterale idrofobico
                res1_group = f"sc.{res1_name.lower()}"
                res2_group = f"sc.{res2_name.lower()}"

                persistence_percentage = (persistence / frame_count) * 100

                # Scrivi i dati nel formato richiesto
                f.write(
                    f"SYSTEM\t{res1_id}\t{res1_name}\t{res1_group}\tSYSTEM\t{res2_id}\t{res2_name}\t{res2_group}\t{persistence_percentage:.2f}\n")

        # Salva la matrice di persistenza in un file .dat
        np.savetxt(persistence_matrix_file, persistence_matrix, fmt='%.1f')

        # Messaggio di conferma di salvataggio dei file
        self.label_result.configure(text="Files saved successfully!", text_color="white")
"""

# Avvio della GUI
if __name__ == "__main__":
    ctk.set_appearance_mode("dark")  # Modalità scura
    ctk.set_default_color_theme("blue")  # Imposta un tema di default

    root = ctk.CTk()  # Crea la finestra principale con customtkinter
    app = ContactCalculatorHc(root, "topology.pdb", "trajectory.dcd", "output_directory")
    root.mainloop()
