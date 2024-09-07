import customtkinter as ctk

class PSNEvaluator:
    def __init__(self, root):
        self.root = root
        self.root.title("PSN evaluator")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Crea un Frame esterno con la cornice
        self.outer_frame = ctk.CTkFrame(root, fg_color='lightgreen', border_color="black", border_width=2)
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Crea un Frame interno per contenere i widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color='lightgreen')
        self.frame.pack(expand=True, padx=10, pady=10)

        # Selezione della soglia
        self.label_threshold = ctk.CTkLabel(self.frame, text="Select the threshold:", font=("Helvetica", 50), text_color="black")
        self.label_threshold.grid(row=0, column=0, padx=10, pady=10)

        self.threshold_var = ["0.5", "1.0", "1.5", "2.0"]
        self.selected_threshold_var = ctk.StringVar()
        self.combobox_threshold = ctk.CTkComboBox(self.frame, variable=self.selected_threshold_var, values=self.threshold_var, width=150, height=40, font=("Helvetica", 35))
        self.combobox_threshold.grid(row=0, column=1, padx=10, pady=10)
        self.combobox_threshold.set(self.threshold_var[1])  # Valore di default

        # Selezione delle metriche aggiuntive
        self.label_metrics = ctk.CTkLabel(self.frame, text="Select the additional metrics:", font=("Helvetica", 50), text_color="black")
        self.label_metrics.grid(row=1, column=0, padx=10, pady=10)

        self.metrics_var = ["Metric1", "Metric2", "Metric3"]
        self.selected_metrics_var = ctk.StringVar()
        self.combobox_metrics = ctk.CTkComboBox(self.frame, variable=self.selected_metrics_var, values=self.metrics_var, width=150, height=40, font=("Helvetica", 35))
        self.combobox_metrics.grid(row=1, column=1, padx=10, pady=10)
        self.combobox_metrics.set(self.metrics_var[0])  # Valore di default

        # Selezione dei formati di output
        self.label_output = ctk.CTkLabel(self.frame, text="Select the output formats:", font=("Helvetica", 50), text_color="black")
        self.label_output.grid(row=2, column=0, padx=10, pady=10)

        self.var_gml = ctk.BooleanVar()
        self.var_png = ctk.BooleanVar()
        self.var_pml = ctk.BooleanVar()
        self.var_vmd = ctk.BooleanVar()
        self.var_pdb = ctk.BooleanVar()

        self.check_gml = ctk.CTkCheckBox(self.frame, text=".gml", variable=self.var_gml, font=("Helvetica", 35), text_color="black")
        self.check_png = ctk.CTkCheckBox(self.frame, text=".png", variable=self.var_png, font=("Helvetica", 35), text_color="black")
        self.check_pml = ctk.CTkCheckBox(self.frame, text="PyMOL script (.pml)", variable=self.var_pml, font=("Helvetica", 35), text_color="black")
        self.check_vmd = ctk.CTkCheckBox(self.frame, text="VMD script (.vmd)", variable=self.var_vmd, font=("Helvetica", 35), text_color="black")
        self.check_pdb = ctk.CTkCheckBox(self.frame, text=".pdb", variable=self.var_pdb, font=("Helvetica", 35), text_color="black")

        self.check_gml.grid(row=3, column=0, sticky='w', padx=10)
        self.check_png.grid(row=3, column=1, sticky='w', padx=10)
        self.check_pml.grid(row=4, column=0, sticky='w', padx=10)
        self.check_vmd.grid(row=4, column=1, sticky='w', padx=10)
        self.check_pdb.grid(row=5, column=0, sticky='w', padx=10)

        # Bottone per avviare il processo
        self.button_process = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 35), command=self.process_files, width=150, height=40, fg_color="gray30")
        self.button_process.grid(row=6, column=0, columnspan=2, pady=20)

        # Etichetta per il risultato
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 35), text_color="black")
        self.label_result.grid(row=7, column=0, columnspan=2, pady=10)

    def process_files(self):
        threshold = self.selected_threshold_var.get()
        metrics = self.selected_metrics_var.get()
        if not threshold or not metrics:
            self.label_result.config(text="Error: Select a threshold and the additional metrics.", text_color='red')
            return

if __name__ == "__main__":
    root = ctk.CTk()
    app = PSNEvaluator(root)
    root.mainloop()
