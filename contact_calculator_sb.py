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
        self.root.configure(fg_color='#dfff00')  # Colore di sfondo della finestra

        # Crea un frame esterno con bordi e colore di sfondo personalizzato
        self.outer_frame = ctk.CTkFrame(root, fg_color='#dfff00', corner_radius=10)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Crea un frame interno all'interno del frame esterno
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color='#dfff00')
        self.frame.pack(expand=True, padx=10, pady=10)

        # Sezione per inserire la distanza di cutoff per il calcolo dei contatti
        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("Helvetica", 50), fg_color='#dfff00', text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10)

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("Helvetica", 35))
        self.entry_cutOff.insert(0, "4.5")
        self.entry_cutOff.grid(row=2, column=0, padx=10, pady=10)

        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 35), command=self.confirm_cutoff, fg_color="gray30")
        self.button_cutOff.grid(row=2, column=1, padx=10, pady=10)

        # Bottone per avviare il processo di calcolo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 35), command=self.preprocess_file, fg_color="gray30")
        self.button_process.grid(row=3, column=0, columnspan=3, pady=20)

        # Label per visualizzare il risultato dell'operazione (es. completamento o errore)
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 35), fg_color='#dfff00', text_color="black")
        self.label_result.grid(row=4, column=0, columnspan=3, pady=10)

    # Funzione per confermare la distanza di cut-off inserita
    def confirm_cutoff(self):
        self.entry_cutOff.get()

    # Funzione principale per avviare il calcolo, verificando che tutti i dati siano corretti
    def preprocess_file(self):
        topology_file = self.topology_file
        trajectory_file = self.trajectory_file
        cutoff_distance = float(self.entry_cutOff.get())
        output_directory = self.output_directory

        if topology_file and trajectory_file and cutoff_distance:
            self.process_file(topology_file, trajectory_file, cutoff_distance, output_directory)  # Avvia l'elaborazione dei file
            self.label_result.config(text="Processing completed successfully!", text_color="green")  # Mostra un messaggio di successo
        else:
            self.label_result.config(text="Please fill all fields correctly.", text_color="red")  # Messaggio di errore se i campi non sono completi

    # Funzione che esegue effettivamente il calcolo dei ponti salini
    def process_file(self, topology_file, trajectory_file, cutoff_distance, output_directory):
        u = mda.Universe(topology_file, trajectory_file)  # Carica la topologia e la traiettoria usando MDAnalysis

        # Ottieni il nome del file di topologia senza l'estensione
        file_prefix = os.path.splitext(os.path.basename(topology_file))[0]

        # Funzione per trovare i ponti salini utilizzando KD-tree e le distanze tra residui acidi e basici
        def find_salt_bridges(universe, cutoff=3.2):
            acidic_residues = universe.select_atoms(
                "resname ASP GLU and name OE1 OE2 OD1 OD2")  # Seleziona residui acidi
            basic_residues = universe.select_atoms("resname ARG LYS and name NH1 NH2 NZ")  # Seleziona residui basici
            n_residues = len(universe.residues)

            contacts = []  # Lista per memorizzare i contatti tra residui
            adjacency_matrix = np.zeros((n_residues, n_residues))  # Matrice di persistenza inizializzata a 0
            frame_count = 0

            # Loop per ogni frame della traiettoria
            for ts in universe.trajectory:
                frame_count += 1  # Conta i frame
                # Crea KD-tree per una ricerca spaziale efficiente
                kd_tree_acidic = cKDTree(acidic_residues.positions)
                kd_tree_basic = cKDTree(basic_residues.positions)

                # Trova coppie di residui entro il raggio di cutoff
                pairs = kd_tree_acidic.query_ball_tree(kd_tree_basic, r=cutoff)

                # Aggiorna i contatti e la matrice di adiacenza
                for i, connected in enumerate(pairs):
                    for j in connected:
                        acid_res = acidic_residues[i].residue
                        basic_res = basic_residues[j].residue
                        contacts.append((acid_res, basic_res))
                        adjacency_matrix[acid_res.resid, basic_res.resid] += 1  # Aggiungi una persistenza
                        adjacency_matrix[basic_res.resid, acid_res.resid] += 1

            # Converti la matrice di adiacenza in persistenza (percentuale di frame in cui il contatto esiste)
            persistence_matrix = (adjacency_matrix / frame_count) * 100

            return contacts, persistence_matrix  # Ritorna i contatti e la matrice di persistenza

        # Chiama la funzione per trovare i ponti salini
        contacts, persistence_matrix = find_salt_bridges(u, cutoff=cutoff_distance)

        # Salva i risultati: i contatti vengono salvati in un file di testo e la matrice di persistenza in un file .dat
        contacts_file = os.path.join(output_directory, f"{file_prefix}_contacts_SB.txt")
        persistence_matrix_file = os.path.join(output_directory, f"{file_prefix}_matrix_SB.dat")

        with open(contacts_file, 'w') as f:
            f.write("Residuo1\tResiduo2\tPersistenza\n")
            for i in range(len(persistence_matrix)):
                for j in range(i + 1, len(persistence_matrix)):
                    if persistence_matrix[i, j] > 0:  # Solo contatti con persistenza
                        acid_res, basic_res = contacts[i]
                        f.write(f"{acid_res.resid}\t{basic_res.resid}\t{persistence_matrix[i, j]:.1f}\n")

        np.savetxt(persistence_matrix_file, persistence_matrix, fmt='%.1f')

        self.label_result.configure(text="Files saved successfully!", text_color="green")  # Messaggio di conferma di salvataggio dei file


# Avvio della GUI
if __name__ == "__main__":
    ctk.set_appearance_mode("dark")  # Imposta la modalità chiara (opzionale)
    ctk.set_default_color_theme("green")  # Imposta il tema di default (opzionale)

    root = ctk.CTk()  # Usa CTk per creare la finestra principale
    app = ContactCalculatorSb(root, "topology.pdb", "trajectory.dcd", "output_directory")
    root.mainloop()
