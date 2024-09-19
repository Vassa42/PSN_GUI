import customtkinter as ctk
import numpy as np
import MDAnalysis as mda
from scipy.spatial import cKDTree
import os

class ContactCalculatorSb:
    def __init__(self, root, topology_file, trajectory_file, output_directory):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.output_directory = output_directory
        self.root = root
        self.root.title("Contact Calculator: Salt Bridges")  # Titolo della finestra GUI

        # Crea un frame esterno con bordi e colore di sfondo personalizzato
        self.outer_frame = ctk.CTkFrame(root, fg_color='SteelBlue2', corner_radius=10)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Crea un frame interno all'interno del frame esterno
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color='SteelBlue2')
        self.frame.pack(expand=True, padx=10, pady=10)

        # Sezione per inserire la distanza di cutoff per il calcolo dei contatti
        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("Helvetica", 40), fg_color='SteelBlue2', text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("Helvetica", 30))
        self.entry_cutOff.insert(0, "4.5")
        self.entry_cutOff.grid(row=1, column=1, padx=10, pady=10)

        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30), command=self.confirm_cutoff, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=1, column=2, padx=10, pady=10)

        self.label_charged = ctk.CTkLabel(self.frame, text="Select the charged groups' file:", font=("Helvetica", 40), text_color="black")
        self.label_charged.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        self.entry_charged = ctk.CTkEntry(self.frame, width=500, justify="center", font=("Helvetica", 25))
        self.entry_charged.insert(0, "/home/vassa42/Documenti/Tirocinio-Tesi/PSN_GUI/charged_groups.ini")
        self.entry_charged.grid(row=2, column=1, padx=10, pady=10)

        self.button_charged = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 30), command=self.browse_charged_groups_file, width=200, height=40, fg_color="gray30")
        self.button_charged.grid(row=2, column=2, padx=10, pady=10)

        # Bottone per avviare il processo di calcolo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), fg_color="gray30", width=200)
        self.button_process.grid(row=3, column=0, columnspan=3, pady=20)

        # Label per visualizzare il risultato dell'operazione (es. completamento o errore)
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40), fg_color='SteelBlue2', text_color="black")
        self.label_result.grid(row=4, column=0, columnspan=3, pady=10)

    # Funzione per confermare la distanza di cut-off inserita
    def confirm_cutoff(self):
        self.entry_cutOff.get()

    def browse_charged_groups_file(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the charged groups' file", filetypes=[("INI files", "*.ini"), ("All files", "*.*")])
        if filename:
            self.entry_charged.delete(0, ctk.END)
            self.entry_charged.insert(0, filename)

""" Funzione principale per avviare il calcolo
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
        u = mda.Universe(topology_file, trajectory_file)

        file_prefix = os.path.splitext(os.path.basename(topology_file))[0]

        def find_salt_bridges(universe, cutoff=3.2):
            acidic_residues = universe.select_atoms("resname ASP GLU and name OE1 OE2 OD1 OD2")
            basic_residues = universe.select_atoms("resname ARG LYS and name NH1 NH2 NZ")
            n_residues = len(universe.residues)

            contacts = []
            adjacency_matrix = np.zeros((n_residues, n_residues))
            frame_count = 0

            # Nuovo set per evitare doppi conteggi di contatti nello stesso frame
            contacts_set = set()

            for ts in universe.trajectory:
                frame_count += 1
                kd_tree_acidic = cKDTree(acidic_residues.positions)
                kd_tree_basic = cKDTree(basic_residues.positions)
                pairs = kd_tree_acidic.query_ball_tree(kd_tree_basic, r=cutoff)

                for i, connected in enumerate(pairs):
                    for j in connected:
                        acid_res = acidic_residues[i].residue
                        basic_res = basic_residues[j].residue

                        # Creiamo una chiave unica per evitare doppi contatti nello stesso frame
                        contact_key = tuple(sorted((acid_res.resid, basic_res.resid)))
                        if contact_key not in contacts_set:
                            contacts.append((acid_res, basic_res))
                            adjacency_matrix[acid_res.resid, basic_res.resid] += 1
                            adjacency_matrix[basic_res.resid, acid_res.resid] += 1
                            contacts_set.add(contact_key)

                contacts_set.clear()  # Reset dei contatti per il prossimo frame

            # Persistenza come percentuale di durata del legame sul totale dei frame
            persistence_matrix = (adjacency_matrix / frame_count) * 100
            return contacts, persistence_matrix

        contacts, persistence_matrix = find_salt_bridges(u, cutoff=cutoff_distance)

        # Salva i risultati dei contatti nel formato richiesto
        contacts_file = os.path.join(output_directory, f"{file_prefix}_contacts_SB.txt")

        with open(contacts_file, 'w') as f:
            f.write(
                "#CHAIN\tRES1_ID\tRES1_NAME\tRES1_CHARGEDGROUP\tRES2_ID\tRES2_NAME\tRES2_CHARGEGROUP\tOCCURRENCE_PERC\n")

            for i in range(len(persistence_matrix)):
                for j in range(i + 1, len(persistence_matrix)):
                    if persistence_matrix[i, j] > 0:
                        acid_res, basic_res = contacts[i]

                        chain1 = acid_res.segid
                        chain2 = basic_res.segid
                        res1_id = acid_res.resid
                        res1_name = acid_res.resname
                        res2_id = basic_res.resid
                        res2_name = basic_res.resname

                        # Gruppi caricati (puoi personalizzare ulteriormente)
                        if res1_name == "ASP":
                            res1_group = "sc.asp.COOn"
                        elif res1_name == "GLU":
                            res1_group = "sc.glu.COOn"
                        if res2_name == "LYS":
                            res2_group = "sc.lys.NZp"
                        elif res2_name == "ARG":
                            res2_group = "sc.arg.NH1p"

                        # Scrivi i dati nel formato richiesto
                        f.write(
                            f"SYSTEM\t{res1_id}\t{res1_name}\t{res1_group}\tSYSTEM\t{res2_id}\t{res2_name}\t{res2_group}\t{persistence_matrix[i, j]:.2f}\n")

        # Salva la matrice di persistenza in un file .dat
        persistence_matrix_file = os.path.join(output_directory, f"{file_prefix}_matrix_SB.dat")
        np.savetxt(persistence_matrix_file, persistence_matrix, fmt='%.1f')

        self.label_result.configure(text="Files saved successfully!", text_color="white")
"""
# Avvio della GUI
if __name__ == "__main__":
    ctk.set_appearance_mode("dark")
    ctk.set_default_color_theme("blue")

    root = ctk.CTk()
    app = ContactCalculatorSb(root, "topology.pdb", "trajectory.dcd", "output_directory")
    root.mainloop()
