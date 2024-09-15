import customtkinter as ctk
import numpy as np
import MDAnalysis as mda
import os
from scipy.spatial import cKDTree


class ContactCalculatorHb:
    def __init__(self, root, topology_file, trajectory_file, output_directory):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.output_directory = output_directory
        self.root = root
        self.root.title("Contact Calculator Hydrogen Bonds")

        # Imposta tema scuro di default
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")  # Tema colore verde

        # Frame esterno con bordi e padding
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="DodgerBlue2")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="DodgerBlue2")
        self.frame.pack(expand=True, padx=10, pady=10)

        # Label e campo per l'inserimento della distanza di cut-off
        self.label_cutOff = ctk.CTkLabel(self.frame, text="Insert the cut-off distance (A):", font=("helvetica", 50),
                                         text_color="black")
        self.label_cutOff.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.entry_cutOff = ctk.CTkEntry(self.frame, width=200, justify="center", font=("helvetica", 40))
        self.entry_cutOff.insert(0, "3.5")
        self.entry_cutOff.grid(row=1, column=1, padx=10, pady=10)

        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40),
                                           command=self.confirm_cutoff, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=1, column=2, padx=10, pady=10)

        # Label e campo per l'inserimento dell'angolo del legame
        self.label_angle = ctk.CTkLabel(self.frame, text="Insert the bond's angle (D):",
                                        font=("helvetica", 50), text_color="black")
        self.label_angle.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        self.entry_angle = ctk.CTkEntry(self.frame, width=200, justify="center", font=("helvetica", 40))
        self.entry_angle.insert(0, "120.0")
        self.entry_angle.grid(row=2, column=1, padx=10, pady=10)

        self.button_angle = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40),
                                          command=self.confirm_angle, fg_color="gray30", width=200)
        self.button_angle.grid(row=2, column=2, padx=10, pady=10)

        # Label e campo per la selezione della classe di legami idrogeno
        self.label_class = ctk.CTkLabel(self.frame, text="Select the class to analyze:",
                                        font=("Helvetica", 50), text_color="black")
        self.label_class.grid(row=3, column=0, padx=10, pady=10, sticky='w')

        self.classes = ["All", "Mc-Mc", "Mc-Sc", "Sc-Sc"]
        self.selected_class = ctk.StringVar()

        # Menu a tendina con il font personalizzato
        self.combobox_class = ctk.CTkComboBox(self.frame,
                                              values=self.classes,
                                              variable=self.selected_class,
                                              state='readonly',
                                              width=200,
                                              font=("helvetica", 40),
                                              dropdown_font=(
                                              "helvetica", 40))  # Imposta il font per i valori nel menu a tendina
        self.combobox_class.grid(row=3, column=1, padx=10, pady=10)
        self.combobox_class.set(self.classes[0])

        # Bottone per avviare il processo di calcolo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 40),
                                            command=self.process_files, fg_color="gray30", width=200)
        self.button_process.grid(row=4, column=0, columnspan=3, pady=20)

        # Label per visualizzare il risultato o eventuali messaggi di errore
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40))
        self.label_result.grid(row=5, column=0, columnspan=3, pady=10)

    def confirm_cutoff(self):
        self.entry_cutOff.get()

    def confirm_angle(self):
        self.entry_angle.get()

    def process_files(self):
        topology_file = self.topology_file
        trajectory_file = self.trajectory_file
        cutoff = float(self.entry_cutOff.get())
        hb_class = self.selected_class.get()
        output_directory = self.output_directory

        if not output_directory:
            self.label_result.configure(text="Please select an output directory.", text_color="red")
            return

        # Estrai il nome base del file di topologia senza estensione
        base_filename = os.path.splitext(os.path.basename(topology_file))[0]

        # Carica la topologia e la traiettoria usando MDAnalysis
        try:
            u = mda.Universe(topology_file, trajectory_file)
        except Exception as e:
            self.label_result.configure(text=f"Error loading files: {e}")
            return

        # Seleziona gli atomi e i criteri per i legami idrogeno in base alla classe selezionata
        if hb_class == "Mc-Mc":
            donors = u.select_atoms("protein and name N H")
            acceptors = u.select_atoms("protein and name O")
        elif hb_class == "Mc-Sc":
            donors = u.select_atoms("protein and name N H")
            acceptors = u.select_atoms("protein and not name O")
        elif hb_class == "Sc-Sc":
            donors = u.select_atoms("protein and not name N H")
            acceptors = u.select_atoms("protein and not name O")
        else:  # All
            donors = u.select_atoms("protein and name N H")
            acceptors = u.select_atoms("protein and name O")

        if len(donors) == 0 or len(acceptors) == 0:
            self.label_result.configure(text="No donors or acceptors found based on the selection.")
            return

        # Inizializza la matrice di persistenza (continua)
        n_donors = len(donors)
        n_acceptors = len(acceptors)
        persistence_matrix = np.zeros((n_donors, n_acceptors), dtype=float)

        hb_results = []
        n_frames = len(u.trajectory)

        for ts in u.trajectory:
            donor_positions = donors.positions
            acceptor_positions = acceptors.positions

            # Crea i KD-tree per la ricerca spaziale
            kd_tree_donors = cKDTree(donor_positions)
            kd_tree_acceptors = cKDTree(acceptor_positions)

            pairs = kd_tree_donors.query_ball_tree(kd_tree_acceptors, cutoff)

            for i, acc in enumerate(pairs):
                for j in acc:
                    donor_atom = donors[i]
                    acceptor_atom = acceptors[j]

                    distance = np.linalg.norm(donor_atom.position - acceptor_atom.position)
                    if distance <= cutoff:
                        persistence_matrix[i, j] += 1

                        # Registra i risultati dettagliati
                        hb_results.append(
                            f"Frame {ts.frame}: Donor {donor_atom.resname}-{donor_atom.resid}-{donor_atom.name} "
                            f"- Acceptor {acceptor_atom.resname}-{acceptor_atom.resid}-{acceptor_atom.name}, "
                            f"Distance: {distance:.3f} Å")

        if not hb_results:
            self.label_result.configure(text="No hydrogen bonds were found.", text_color="red")
            return

        # Converti la matrice di persistenza in percentuali
        persistence_matrix = (persistence_matrix / n_frames) * 100

        # Salva i risultati in formato .txt e .dat
        output_txt = os.path.join(output_directory, f"{base_filename}_contacts_HB.txt")
        output_dat = os.path.join(output_directory, f"{base_filename}_matrix_HB.dat")

        with open(output_txt, "w") as f_txt:
            f_txt.write(f"Hydrogen Bond Persistence (in % over {n_frames} frames)\n")

            # Funzione per associare residuo e gruppo chimico
            def get_group_info(resname, atom_name):
                if resname in ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"]:
                    return f"sc.{resname.lower()}"
                elif atom_name == "N":
                    return "bb.NH"
                elif atom_name == "O":
                    return "bb.CO"
                else:
                    return f"sc.{resname.lower()}"

            for i in range(n_donors):
                for j in range(n_acceptors):
                    if persistence_matrix[i, j] > 0:
                        donor = donors[i]
                        acceptor = acceptors[j]

                        donor_info = f"A_{donor.resid}{donor.resname}_{get_group_info(donor.resname, donor.name)}"
                        acceptor_info = f"A_{acceptor.resid}{acceptor.resname}_{get_group_info(acceptor.resname, acceptor.name)}"

                        # Scrivi i risultati nel file con il formato richiesto
                        f_txt.write(f"{donor_info}: {acceptor_info}          {persistence_matrix[i, j]:.1f}%\n")

        # Salva la matrice di persistenza in formato .dat
        np.savetxt(output_dat, persistence_matrix, fmt="%.1f")

        self.label_result.configure(text="Files saved successfully!", text_color="white")


# Avvio della GUI
if __name__ == "__main__":
    root = ctk.CTk()
    app = ContactCalculatorHb(root, "topology.psf", "trajectory.dcd", "output_directory")
    root.mainloop()