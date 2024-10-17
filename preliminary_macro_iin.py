import customtkinter as ctk

class PreliminaryMacroIIN:
    def __init__(self, root):
        self.root = root
        self.root.title("Preliminary Macro-IIN")  # Imposta il titolo della finestra

        # Imposta il tema e i colori di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="gold2")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per posizionare i widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="gold2")
        self.frame.pack(expand=True, padx=10, pady=10)

        # Label e campo di inserimento per la soglia (threshold)
        self.label_threshold = ctk.CTkLabel(self.frame, text="Insert threshold (%):", font=("Helvetica", 40),
                                            text_color="black")
        self.label_threshold.grid(row=0, column=0, padx=10, pady=10)

        self.entry_threshold = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 30), justify="center")
        self.entry_threshold.grid(row=0, column=1, padx=10, pady=10)
        self.entry_threshold.insert(0, "20")  # Inserisce un valore predefinito nel campo di testo

        # Pulsante per eseguire il processo con la soglia inserita
        self.process_button = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), width=200, height=40,
                                            fg_color="gray30")
        self.process_button.grid(row=1, column=0, padx=10, pady=20)


# Avvio dell'applicazione
if __name__ == "__main__":
    root = ctk.CTk()
    app = PreliminaryMacroIIN(root)
    root.mainloop()
