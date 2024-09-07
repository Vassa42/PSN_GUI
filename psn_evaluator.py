import tkinter as tk
from tkinter import ttk

class PSNEvaluator:
    def __init__(self, root):
        self.root = root
        self.root.title("PSN evaluator")

        # Imposta il colore di sfondo della finestra principale
        self.root.configure(bg='lightgreen')

        # Crea un Frame esterno con la cornice
        self.outer_frame = tk.Frame(root, bg='lightgreen', highlightbackground="black", highlightthickness=2, bd=2, relief='solid')
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Crea un Frame interno per contenere i widget
        self.frame = tk.Frame(self.outer_frame, bg='lightgreen')
        self.frame.pack(expand=True, padx=10, pady=10)

        # Selezione della soglia
        self.label_threshold = tk.Label(self.frame, text="Select the threshold:", font=("Helvetica", 20), bg='lightgreen')
        self.label_threshold.grid(row=0, column=0, padx=10, pady=10)

        self.threshold_var = ["0.5", "1.0", "1.5", "2.0"]
        self.selected_threshold_var = tk.StringVar()
        self.combobox_threshold = ttk.Combobox(self.frame, textvariable=self.selected_threshold_var, values=self.threshold_var, state='readonly', width=10)
        self.combobox_threshold.grid(row=0, column=1, padx=10, pady=10)
        self.combobox_threshold.set(self.threshold_var[1]) # Valore di default

        # Selezione delle metriche aggiuntive
        self.label_metrics = tk.Label(self.frame, text="Select the additional metrics:", font=("Helvetica", 20), bg='lightgreen')
        self.label_metrics.grid(row=1, column=0, padx=10, pady=10)

        self.metrics_var = ["Metric1", "Metric2", "Metric3"]
        self.selected_metrics_var = tk.StringVar()
        self.combobox_metrics = ttk.Combobox(self.frame, textvariable=self.selected_metrics_var, values=self.metrics_var, state='readonly', width=10)
        self.combobox_metrics.grid(row=1, column=1, padx=10, pady=10)
        self.combobox_metrics.set(self.metrics_var[0]) # Valore di default

        # Selezione dei formati di output
        self.label_output = tk.Label(self.frame, text="Select the output formats:", font=("Helvetica", 20), bg='lightgreen')
        self.label_output.grid(row=2, column=0, padx=10, pady=10)

        self.var_gml = tk.BooleanVar()
        self.var_png = tk.BooleanVar()
        self.var_pml = tk.BooleanVar()
        self.var_vmd = tk.BooleanVar()
        self.var_pdb = tk.BooleanVar()

        self.check_gml = tk.Checkbutton(self.frame, text=".gml", variable=self.var_gml, font=("Helvetica", 15), bg='lightgreen')
        self.check_png = tk.Checkbutton(self.frame, text=".png", variable=self.var_png, font=("Helvetica", 15), bg='lightgreen')
        self.check_pml = tk.Checkbutton(self.frame, text="PyMOL script (.pml)", variable=self.var_pml, font=("Helvetica", 15), bg='lightgreen')
        self.check_vmd = tk.Checkbutton(self.frame, text="VMD script (.vmd)", variable=self.var_vmd, font=("Helvetica", 15), bg='lightgreen')
        self.check_pdb = tk.Checkbutton(self.frame, text=".pdb", variable=self.var_pdb, font=("Helvetica", 15), bg='lightgreen')

        self.check_gml.grid(row=3, column=0, sticky='w', padx=10)
        self.check_png.grid(row=3, column=1, sticky='w', padx=10)
        self.check_pml.grid(row=4, column=0, sticky='w', padx=10)
        self.check_vmd.grid(row=4, column=1, sticky='w', padx=10)
        self.check_pdb.grid(row=5, column=0, sticky='w', padx=10)

        # Bottone per avviare il processo
        self.button_process = tk.Button(self.frame, text="Process", font=("Helvetica", 15), command=self.process_files)
        self.button_process.grid(row=6, column=0, columnspan=2, pady=20)

        # Etichetta per il risultato
        self.label_result = tk.Label(self.frame, text="", font=("Helvetica", 15), bg='lightgreen')
        self.label_result.grid(row=7, column=0, columnspan=2, pady=10)

    def process_files(self):
        threshold = self.selected_threshold_var.get()
        metrics = self.selected_metrics_var.get()
        if not threshold or not metrics:
            self.label_result.config(text="Error: Select a threshold and the additional metrics.", fg='red')
            return

if __name__ == "__main__":
    root = tk.Tk()
    app = PSNEvaluator(root)
    root.mainloop()

