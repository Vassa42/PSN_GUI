import customtkinter
import customtkinter as ctk
from tkinter import filedialog

from path_analysis import PathAnalysis
from topology_analysis import TopologyAnalysis


class PSNEvaluator:
    def __init__(self, root):
        self.root = root
        self.root.title("PSN Evaluator")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Crea un Frame esterno con la cornice
        self.outer_frame = ctk.CTkFrame(root, fg_color='dark sea green', border_color="black", border_width=2)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Crea un Frame interno per contenere i widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color='dark sea green')
        self.frame.pack(expand=True, padx=10, pady=10)

        # Selezione del file macro_IIN
        self.label_macro_IIN_file = ctk.CTkLabel(self.frame, text="Select the macro IIN file:", font=("Helvetica", 50), text_color="black")
        self.label_macro_IIN_file.grid(row=0, column=0, padx=10, pady=10, sticky='w')

        self.entry_macro_IIN_file = ctk.CTkEntry(self.frame,font=("Helvetica", 30), width=500)
        self.entry_macro_IIN_file.grid(row=0, column=1, padx=10, pady=10)

        self.button_browse_macro_IIN = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 40), command=self.browse_macro_IIN_file, width=200, fg_color="gray30")
        self.button_browse_macro_IIN.grid(row=0, column=2, padx=10, pady=10)

        # Selezione del file di struttura
        self.label_structure_file = ctk.CTkLabel(self.frame, text="Select the structure file:", font=("Helvetica", 50), text_color="black")
        self.label_structure_file.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.entry_structure_file = ctk.CTkEntry(self.frame, font=("Helvetica", 30), width=500)
        self.entry_structure_file.grid(row=1, column=1, padx=10, pady=10)

        self.button_browse_structure = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 40), command=self.browse_structure_file, width=200, fg_color="gray30")
        self.button_browse_structure.grid(row=1, column=2, padx=10, pady=10)

        # Selezione delle metriche da valutare
        self.label_metrics = ctk.CTkLabel(self.frame, text="Select the metrics to evaluate:", font=("Helvetica", 50), text_color="black")
        self.label_metrics.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        self.var_path_analysis = ctk.BooleanVar()
        self.var_topology_analysis = ctk.BooleanVar()

        self.check_path_analysis = ctk.CTkCheckBox(self.frame, text="Path Analysis", variable=self.var_path_analysis, font=("Helvetica", 40), text_color="black")
        self.check_topology_analysis = ctk.CTkCheckBox(self.frame, text="Topology Analysis", variable=self.var_topology_analysis, font=("Helvetica", 40), text_color="black")

        self.check_path_analysis.grid(row=3, column=0, sticky='w', padx=10)
        self.check_topology_analysis.grid(row=4, column=0, sticky='w', padx=10)

        # Bottone per la selezione
        self.button_process = ctk.CTkButton(self.frame, text="Select", command=self.open_metrics, font=("Helvetica", 40), width=200, height=40, fg_color="gray30")
        self.button_process.grid(row=7, column=0, columnspan=2, pady=20)

        # Etichetta per il risultato
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40), text_color="black")
        self.label_result.grid(row=8, column=0, columnspan=2, pady=10)

    def browse_macro_IIN_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Dat files", "*.dat")])
        if file_path:
            self.entry_macro_IIN_file.delete(0, ctk.END)
            self.entry_macro_IIN_file.insert(0, file_path)

    def browse_structure_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("Pdb files", "*.pdb"), ("Gro files", "*.gro")])
        if file_path:
            self.entry_structure_file.delete(0, ctk.END)
            self.entry_structure_file.insert(0, file_path)

    def open_metrics(self):
        if self.var_path_analysis.get() and self.var_topology_analysis.get():
            self.new_window1 = customtkinter.CTkToplevel(self.root)
            PathAnalysis(self.new_window1)
            self.new_window2 = customtkinter.CTkToplevel(self.root)
            TopologyAnalysis(self.new_window2)
        elif self.var_path_analysis.get():
            self.new_window3 = customtkinter.CTkToplevel(self.root)
            PathAnalysis(self.new_window3)
        elif self.var_topology_analysis.get():
            self.new_window4 = customtkinter.CTkToplevel(self.root)
            TopologyAnalysis(self.new_window4)
"""
    def process_files(self):
        macro_IIN_file = self.entry_macro_IIN_file.get()
        structure_file = self.entry_structure_file.get()
        graph_analysis = self.var_graph_analysis.get()
        path_analysis = self.var_path_analysis.get()
        centrality_analysis = self.var_centrality_analysis.get()

        if not macro_IIN_file or not structure_file:
            self.label_result.config(text="Error: Select both macro IIN file and structure file.", text_color='red')
            return

        selected_metrics = []
        if graph_analysis:
            selected_metrics.append("Graph Analysis")
        if path_analysis:
            selected_metrics.append("Path Analysis")
        if centrality_analysis:
            selected_metrics.append("Centrality Analysis")

        if not selected_metrics:
            self.label_result.config(text="Error: Select at least one metric.", text_color='red')
            return

        # Messaggio di conferma, per ora non implementa la logica del processo
        self.label_result.config(text=f"Processing with macro IIN file: {macro_IIN_file}\nFile: {structure_file}\nMetrics: {', '.join(selected_metrics)}", text_color='white')
"""

if __name__ == "__main__":
    root = ctk.CTk()
    app = PSNEvaluator(root)
    root.mainloop()
