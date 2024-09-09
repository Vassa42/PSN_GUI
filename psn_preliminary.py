import customtkinter as ctk
from customtkinter import filedialog
import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt

class PSNPreliminary:
    def __init__(self, root):
        self.root = root
        self.root.title("PSN Preliminary")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="orange")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per il posizionamento dei widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="orange")
        self.frame.pack(expand=True, padx=10, pady=10)

        # Valori Min, Max e Step per il filtraggio
        self.label_min = ctk.CTkLabel(self.frame, text="Minimum value:", font=("Helvetica", 50), text_color="black")
        self.label_min.grid(row=0, column=0, padx=10, pady=10)

        self.entry_min = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 40))
        self.entry_min.grid(row=0, column=1, padx=10, pady=10)
        self.entry_min.insert(0, "0.0")

        self.button_min = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40), command=self.update_min, width=200, height=40, fg_color="gray30")
        self.button_min.grid(row=0, column=2, padx=10, pady=10)

        self.label_max = ctk.CTkLabel(self.frame, text="Maximum value:", font=("Helvetica", 50), text_color="black")
        self.label_max.grid(row=1, column=0, padx=10, pady=10)

        self.entry_max = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 40))
        self.entry_max.grid(row=1, column=1, padx=10, pady=10)
        self.entry_max.insert(0, "100.0")

        self.button_max = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40), command=self.update_max, width=200, height=40, fg_color="gray30")
        self.button_max.grid(row=1, column=2, padx=10, pady=10)

        self.label_step = ctk.CTkLabel(self.frame, text="Step value:", font=("Helvetica", 50), text_color="black")
        self.label_step.grid(row=2, column=0, padx=10, pady=10)

        self.entry_step = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 40))
        self.entry_step.grid(row=2, column=1, padx=10, pady=10)
        self.entry_step.insert(0, "5.0")

        self.button_step = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40), command=self.update_step, width=200, height=40, fg_color="gray30")
        self.button_step.grid(row=2, column=2, padx=10, pady=10)

        # Checkbox per soglia manuale
        self.threshold_var = ctk.IntVar()
        self.checkbox_threshold = ctk.CTkCheckBox(self.frame, text="Use manual threshold", font=("Helvetica", 40), variable=self.threshold_var, text_color="black")
        self.checkbox_threshold.grid(row=3, column=0, padx=10, pady=10)

        self.label_threshold = ctk.CTkLabel(self.frame, text="Threshold:", font=("Helvetica", 40), text_color="black")
        self.label_threshold.grid(row=3, column=1, padx=10, pady=10)

        self.entry_threshold = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 40))
        self.entry_threshold.grid(row=3, column=2, padx=10, pady=10)
        self.entry_threshold.insert(0, "0.5")

        # Inserimento file matrici di persistenza
        self.label_persistence_SB = ctk.CTkLabel(self.frame, text="Persistence matrix SB:", font=("Helvetica", 50), text_color="black")
        self.label_persistence_SB.grid(row=4, column=0, padx=10, pady=10)

        self.entry_persistence_SB = ctk.CTkEntry(self.frame, width=200, font=("Helvetica", 40))
        self.entry_persistence_SB.grid(row=4, column=1, padx=10, pady=10)

        self.button_persistence_SB = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 40), command=self.browse_persistence_file_SB, width=200, height=40, fg_color="gray30")
        self.button_persistence_SB.grid(row=4, column=2, padx=10, pady=10)

        self.label_persistence_HC = ctk.CTkLabel(self.frame, text="Persistence matrix HC:", font=("Helvetica", 50), text_color="black")
        self.label_persistence_HC.grid(row=5, column=0, padx=10, pady=10)

        self.entry_persistence_HC = ctk.CTkEntry(self.frame, width=200, font=("Helvetica", 40))
        self.entry_persistence_HC.grid(row=5, column=1, padx=10, pady=10)

        self.button_persistence_HC = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 40), command=self.browse_persistence_file_HC, width=200, height=40, fg_color="gray30")
        self.button_persistence_HC.grid(row=5, column=2, padx=10, pady=10)

        self.label_persistence_HB = ctk.CTkLabel(self.frame, text="Persistence matrix HB:", font=("Helvetica", 50), text_color="black")
        self.label_persistence_HB.grid(row=6, column=0, padx=10, pady=10)

        self.entry_persistence_HB = ctk.CTkEntry(self.frame, width=200, font=("Helvetica", 40))
        self.entry_persistence_HB.grid(row=6, column=1, padx=10, pady=10)

        self.button_persistence_HB = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 40), command=self.browse_persistence_file_HB, width=200, height=40, fg_color="gray30")
        self.button_persistence_HB.grid(row=6, column=2, padx=10, pady=10)

        # Checkbox per generare Macro IIN e Plot
        self.macro_iin_var = ctk.IntVar()
        self.checkbox_macro_iin = ctk.CTkCheckBox(self.frame, text="Generate Macro IIN", font=("Helvetica", 40),
                                                  variable=self.macro_iin_var, text_color="black")
        self.checkbox_macro_iin.grid(row=7, column=0, padx=10, pady=10)

        self.plot_var = ctk.IntVar()
        self.checkbox_plot = ctk.CTkCheckBox(self.frame, text="Generate Plot", font=("Helvetica", 40),
                                             variable=self.plot_var, text_color="black")
        self.checkbox_plot.grid(row=8, column=0, padx=10, pady=10)

        # Pulsante per avviare il processo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 40), command=self.process_files, width=200, height=40, fg_color="gray30")
        self.button_process.grid(row=9, column=0, columnspan=3, pady=20)

        # Label per il risultato
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40))
        self.label_result.grid(row=10, column=0, columnspan=3, pady=10)

        # Valori di default
        self.min_value = 0.0
        self.max_value = 100.0
        self.step_value = 5.0

    def update_min(self):
        try:
            self.min_value = float(self.entry_min.get())
        except ValueError:
            self.label_result.configure(text="Error: Invalid min value.", text_color="red")

    def update_max(self):
        try:
            self.max_value = float(self.entry_max.get())
        except ValueError:
            self.label_result.configure(text="Error: Invalid max value.", text_color="red")

    def update_step(self):
        try:
            self.step_value = float(self.entry_step.get())
        except ValueError:
            self.label_result.configure(text="Error: Invalid step value.", text_color="red")

    # Funzioni di caricamento file
    def browse_persistence_file_SB(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the persistence matrix SB",
                                                  filetypes=[("DAT files", "*.dat"), ("All files", "*.*")])
        if filename:
            self.entry_persistence_SB.delete(0, ctk.END)
            self.entry_persistence_SB.insert(0, filename)

    def browse_persistence_file_HC(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the persistence matrix HC",
                                                  filetypes=[("DAT files", "*.dat"), ("All files", "*.*")])
        if filename:
            self.entry_persistence_HC.delete(0, ctk.END)
            self.entry_persistence_HC.insert(0, filename)

    def browse_persistence_file_HB(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the persistence matrix HB",
                                                  filetypes=[("DAT files", "*.dat"), ("All files", "*.*")])
        if filename:
            self.entry_persistence_HB.delete(0, ctk.END)
            self.entry_persistence_HB.insert(0, filename)

    # Funzione per processare i file
    def process_files(self):
        try:
            min_value = self.min_value
            max_value = self.max_value
            step_value = self.step_value
        except ValueError:
            self.label_result.configure(text="Error: Invalid min, max, or step values.", text_color="red")
            return

        persistence_file_SB = self.entry_persistence_SB.get()
        persistence_file_HC = self.entry_persistence_HC.get()
        persistence_file_HB = self.entry_persistence_HB.get()

        if not persistence_file_SB or not persistence_file_HC or not persistence_file_HB:
            self.label_result.configure(text="Error: Please select all the input files.", text_color="red")
            return

        # Load persistence matrices
        try:
            persistence_matrix_SB = self.load_persistence_matrix(persistence_file_SB)
            persistence_matrix_HC = self.load_persistence_matrix(persistence_file_HC)
            persistence_matrix_HB = self.load_persistence_matrix(persistence_file_HB)

            # Resize persistence matrices to match the largest one
            persistence_matrix_SB, persistence_matrix_HC, persistence_matrix_HB = self.resize_matrices(
                persistence_matrix_SB, persistence_matrix_HC, persistence_matrix_HB
            )

            # Ensure all matrices are 2D
            if not (
                    persistence_matrix_SB.ndim == 2 and persistence_matrix_HC.ndim == 2 and persistence_matrix_HB.ndim == 2):
                raise ValueError("All persistence matrices must be 2D.")
        except Exception as e:
            self.label_result.configure(text=f"Error loading persistence matrices: {e}", text_color="red")
            return

        # Select output folder
        output_folder = filedialog.askdirectory(parent=self.root, title="Select output folder")
        if not output_folder:
            self.label_result.configure(text="Error: No output folder selected.", text_color="red")
            return

        base_name = os.path.basename(persistence_file_SB)
        base_name = base_name.split("matrix")[0].rstrip("_")

        try:
            if self.threshold_var.get():
                threshold = float(self.entry_threshold.get())
            else:
                # Find optimal threshold if not provided manually
                threshold = self.find_optimal_threshold(persistence_matrix_SB, persistence_matrix_HC,
                                                        persistence_matrix_HB)
                print(f"Optimal threshold found: {threshold}")

            adj_matrix_SB = (persistence_matrix_SB > threshold).astype(int)
            adj_matrix_HC = (persistence_matrix_HC > threshold).astype(int)
            adj_matrix_HB = (persistence_matrix_HB > threshold).astype(int)

            # Apply boolean OR
            combined_matrix = np.logical_or.reduce([adj_matrix_SB, adj_matrix_HC, adj_matrix_HB])

            # Controlliamo lo stato delle checkbox e chiamiamo i metodi corrispondenti
            if self.macro_iin_var.get():
                self.create_macro_iin(combined_matrix, base_name, output_folder)

            if self.plot_var.get():
                combined_graph = nx.from_numpy_array(combined_matrix)
                self.create_plot(combined_graph, base_name, output_folder)

            self.label_result.configure(text="Processing complete.", text_color="white")

        except Exception as e:
            self.label_result.configure(text=f"Error: {e}", text_color="red")

    # Funzioni ausiliarie per la gestione delle matrici
    def resize_matrices(self, matrix_SB, matrix_HC, matrix_HB):
        max_rows = max(matrix_SB.shape[0], matrix_HC.shape[0], matrix_HB.shape[0])
        max_cols = max(matrix_SB.shape[1], matrix_HC.shape[1], matrix_HB.shape[1])

        # Function to resize a matrix by adding zeros
        def resize_matrix(matrix, target_rows, target_cols):
            rows, cols = matrix.shape
            resized_matrix = np.zeros((target_rows, target_cols))
            resized_matrix[:rows, :cols] = matrix
            return resized_matrix

        resized_SB = resize_matrix(matrix_SB, max_rows, max_cols)
        resized_HC = resize_matrix(matrix_HC, max_rows, max_cols)
        resized_HB = resize_matrix(matrix_HB, max_rows, max_cols)

        return resized_SB, resized_HC, resized_HB

    def find_optimal_threshold(self, matrix_SB, matrix_HC, matrix_HB):
        # Placeholder logic for finding optimal threshold
        combined_matrix = np.logical_or.reduce([matrix_SB, matrix_HC, matrix_HB])
        # Calculate the average value as a basic heuristic for optimal threshold
        optimal_threshold = np.mean(combined_matrix)
        return optimal_threshold

    def load_persistence_matrix(self, file_path):
        return np.loadtxt(file_path)

    def generate_psn_dat(self, matrix, file_path):
        """Save the adjacency matrix to a .dat file."""
        np.savetxt(file_path, matrix, fmt='%d')

    def plot_graph(self, graph, file_path):
        pos = nx.spring_layout(graph)
        plt.figure(figsize=(8, 8))
        nx.draw(graph, pos, node_size=50, with_labels=False)
        plt.savefig(file_path)
        plt.close()

    def plot_clusters_vs_threshold(self, graph, output_folder, base_name):
        # Placeholder for pmin and cluster count computation
        pmin_values = np.linspace(self.min_value, self.max_value, 100)
        cluster_sizes = [len(max(nx.connected_components(graph), key=len)) for _ in pmin_values]

        plt.figure(figsize=(8, 6))
        plt.plot(pmin_values, cluster_sizes, marker='o')
        plt.xlabel('pmin')
        plt.ylabel('Largest Cluster Size')
        plt.title('Threshold vs. Largest Cluster Size')
        plot_output_path = os.path.join(output_folder, f"{base_name}_cluster_vs_pmin.png")
        plt.savefig(plot_output_path)
        plt.close()

    def create_plot(self, graph, base_name, output_folder):
        # Salviamo il grafico in formato .png
        png_output_path = os.path.join(output_folder, f"{base_name}_macroIIN.png")
        self.plot_graph(graph, png_output_path)

        # Plot sigmoid-like graph
        self.plot_clusters_vs_threshold(graph, output_folder, base_name)

    def create_macro_iin(self, combined_matrix, base_name, output_folder):
        # Salviamo la matrice combinata come file .dat
        dat_output_path = os.path.join(output_folder, f"{base_name}_macroIIN.dat")
        self.generate_psn_dat(combined_matrix, dat_output_path)

        # Generiamo e salviamo il grafo in formato .gml
        combined_graph = nx.from_numpy_array(combined_matrix)
        gml_output_path = os.path.join(output_folder, f"{base_name}_macroIIN.gml")
        nx.write_gml(combined_graph, gml_output_path)


if __name__ == "__main__":
    root = ctk.CTk()
    app = PSNPreliminary(root)
    root.mainloop()
