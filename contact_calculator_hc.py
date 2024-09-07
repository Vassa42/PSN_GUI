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
        self.outer_frame = ctk.CTkFrame(root, fg_color="#7fff00", corner_radius=10)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="#7fff00")
        self.frame.pack(expand=True, padx=10, pady=10)

        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("Helvetica", 50), text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10)

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("Helvetica", 35))
        self.entry_cutOff.insert(0, "5.0")
        self.entry_cutOff.grid(row=2, column=0, padx=10, pady=10)

        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 35), command=self.confirm_cutoff, fg_color="gray30")
        self.button_cutOff.grid(row=2, column=1, padx=10, pady=10)

        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 35), command=self.preprocess_file, fg_color="gray30")
        self.button_process.grid(row=3, column=0, columnspan=3, pady=20)

        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 35), text_color="black")
        self.label_result.grid(row=4, column=0, columnspan=3, pady=10)

    def confirm_cutoff(self):
        self.entry_cutOff.get()

    def preprocess_file(self):
        topology_file = self.topology_file
        trajectory_file = self.trajectory_file
        cutoff_distance = float(self.entry_cutOff.get())
        output_directory = self.output_directory

        if topology_file and trajectory_file and cutoff_distance:
            self.process_file(topology_file, trajectory_file, cutoff_distance, output_directory)
            self.label_result.configure(text="Processing completed successfully!", text_color="green")
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

        # Salva i contatti idrofobici in un file di testo
        contacts_file = os.path.join(output_directory, f"{file_prefix}_contacts_HC.txt")
        persistence_matrix_file = os.path.join(output_directory, f"{file_prefix}_matrix_HC.dat")

        with open(contacts_file, 'w') as f:
            f.write("Residuo1\tResiduo2\tPersistenza\n")
            for (i, j), persistence in interaction_dict.items():
                # Solo contatti con persistenza
                res1 = hydrophobic_residues[i]
                res2 = hydrophobic_residues[j]
                res1_info = f"A_{res1.resid}{res1.resname}_sc.{res1.resname.lower()}"
                res2_info = f"A_{res2.resid}{res2.resname}_sc.{res2.resname.lower()}"
                persistence_percentage = (persistence / frame_count) * 100
                f.write(f"{res1_info}: {res2_info}          {persistence_percentage:.1f}\n")

        # Salva la matrice di persistenza in un file .dat
        np.savetxt(persistence_matrix_file, persistence_matrix, fmt='%.1f')

        # Messaggio di conferma di salvataggio dei file
        self.label_result.configure(text="Files saved successfully!", text_color="green")

# Avvio della GUI
if __name__ == "__main__":
    ctk.set_appearance_mode("dark")  # Modalità scura
    ctk.set_default_color_theme("green")  # Imposta un tema di default

    root = ctk.CTk()  # Crea la finestra principale con customtkinter
    app = ContactCalculatorHc(root, "topology.pdb", "trajectory.dcd", "output_directory")
    root.mainloop()
