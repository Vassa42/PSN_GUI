import customtkinter
import customtkinter as ctk

from preliminary_macro_iin import PreliminaryMacroIIN
from preliminary_plot import PreliminaryPlot


class PSNPreliminary:
    def __init__(self, root):
        self.root = root
        self.root.title("PSN Preliminary")  # Titolo della finestra principale

        # Imposta tema e colori
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="orange")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per il posizionamento dei widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="orange")
        self.frame.pack(expand=True, padx=10, pady=10)

        # Sezione per l'inserimento del file di matrice di persistenza per SB (salt bridges)
        self.label_persistence_SB = ctk.CTkLabel(self.frame, text="Persistence matrix SB:", font=("Helvetica", 40), text_color="black")
        self.label_persistence_SB.grid(row=0, column=0, padx=10, pady=10, sticky='w')

        self.entry_persistence_SB = ctk.CTkEntry(self.frame, width=500, font=("Helvetica", 25))
        self.entry_persistence_SB.grid(row=0, column=1, padx=10, pady=10)

        self.button_persistence_SB = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 30), command=self.browse_persistence_file_SB, width=200, height=40, fg_color="gray30")
        self.button_persistence_SB.grid(row=0, column=2, padx=10, pady=10)

        # Sezione per l'inserimento del file di matrice di persistenza per HC (hydrophobic contacts)
        self.label_persistence_HC = ctk.CTkLabel(self.frame, text="Persistence matrix HC:", font=("Helvetica", 40), text_color="black")
        self.label_persistence_HC.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.entry_persistence_HC = ctk.CTkEntry(self.frame, width=500, font=("Helvetica", 25))
        self.entry_persistence_HC.grid(row=1, column=1, padx=10, pady=10)

        self.button_persistence_HC = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 30), command=self.browse_persistence_file_HC, width=200, height=40, fg_color="gray30")
        self.button_persistence_HC.grid(row=1, column=2, padx=10, pady=10)

        # Sezione per l'inserimento del file di matrice di persistenza per HB (hydrogen bonds)
        self.label_persistence_HB = ctk.CTkLabel(self.frame, text="Persistence matrix HB:", font=("Helvetica", 40), text_color="black")
        self.label_persistence_HB.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        self.entry_persistence_HB = ctk.CTkEntry(self.frame, width=500, font=("Helvetica", 25))
        self.entry_persistence_HB.grid(row=2, column=1, padx=10, pady=10)

        self.button_persistence_HB = ctk.CTkButton(self.frame, text="Browse", font=("Helvetica", 30), command=self.browse_persistence_file_HB, width=200, height=40, fg_color="gray30")
        self.button_persistence_HB.grid(row=2, column=2, padx=10, pady=10)

        # Sezione per scegliere l'output da generare (plot e/o macro IIN)
        self.label_choose_output = ctk.CTkLabel(self.frame, text="Output to generate:", font=("Helvetica", 40), text_color="black")
        self.label_choose_output.grid(row=3, column=0, padx=10, pady=10, sticky='w')

        # Checkbox per scegliere di generare un plot
        self.plot_var = ctk.IntVar()
        self.checkbox_plot = ctk.CTkCheckBox(self.frame, text="Generate Plot", font=("Helvetica", 40), variable=self.plot_var, text_color="black")
        self.checkbox_plot.grid(row=4, column=0, padx=10, pady=10, sticky='w')

        # Checkbox per scegliere di generare una Macro IIN (macromolecular interaction network)
        self.macro_iin_var = ctk.IntVar()
        self.checkbox_macro_iin = ctk.CTkCheckBox(self.frame, text="Generate Macro IIN", font=("Helvetica", 40), variable=self.macro_iin_var, text_color="black")
        self.checkbox_macro_iin.grid(row=5, column=0, padx=10, pady=10, sticky='w')

        # Pulsante per avviare la generazione degli output selezionati
        self.button_process = ctk.CTkButton(self.frame, text="Select", font=("Helvetica", 30), command=self.open_preliminaries, width=200, height=40, fg_color="gray30")
        self.button_process.grid(row=6, column=0, columnspan=3, pady=20)

        # Label per visualizzare eventuali risultati o messaggi di errore
        self.label_result = ctk.CTkLabel(self.frame, text="", font=("Helvetica", 40))
        self.label_result.grid(row=7, column=0, columnspan=3, pady=10)

        # Valori di default per eventuali parametri
        self.min_value = 0.0
        self.max_value = 100.0
        self.step_value = 1.0

    # Funzioni per il caricamento dei file di persistenza
    def browse_persistence_file_SB(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the persistence matrix SB", filetypes=[("DAT files", "*.dat"), ("All files", "*.*")])
        if filename:
            self.entry_persistence_SB.delete(0, ctk.END)
            self.entry_persistence_SB.insert(0, filename)

    def browse_persistence_file_HC(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the persistence matrix HC", filetypes=[("DAT files", "*.dat"), ("All files", "*.*")])
        if filename:
            self.entry_persistence_HC.delete(0, ctk.END)
            self.entry_persistence_HC.insert(0, filename)

    def browse_persistence_file_HB(self):
        filename = ctk.filedialog.askopenfilename(parent=self.root, title="Select the persistence matrix HB", filetypes=[("DAT files", "*.dat"), ("All files", "*.*")])
        if filename:
            self.entry_persistence_HB.delete(0, ctk.END)
            self.entry_persistence_HB.insert(0, filename)

    # Funzione per aprire le finestre in base alle scelte (plot e/o macro IIN)
    def open_preliminaries(self):
        if self.plot_var.get() and self.macro_iin_var.get():  # Se entrambe le opzioni sono selezionate
            self.new_window1 = customtkinter.CTkToplevel(self.root)
            PreliminaryPlot(self.new_window1)
            self.new_window2 = customtkinter.CTkToplevel(self.root)
            PreliminaryMacroIIN(self.new_window2)
        elif self.plot_var.get():  # Solo plot
            self.new_window3 = customtkinter.CTkToplevel(self.root)
            PreliminaryPlot(self.new_window3)
        elif self.macro_iin_var.get():  # Solo Macro IIN
            self.new_window4 = customtkinter.CTkToplevel(self.root)
            PreliminaryMacroIIN(self.new_window4)


# Avvio dell'applicazione GUI
if __name__ == "__main__":
    root = ctk.CTk()
    app = PSNPreliminary(root)
    root.mainloop()
