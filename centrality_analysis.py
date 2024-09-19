import customtkinter as ctk


class CentralityAnalysis:
    def __init__(self, root):
        self.root = root
        self.root.title("Centrality Analysis")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="DarkSeaGreen3")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per il posizionamento dei widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="DarkSeaGreen3")
        self.frame.pack(expand=True, padx=10, pady=10)

        self.label_set = ctk.CTkLabel(self.frame, text="Select the set to analyze:",
                                        font=("Helvetica", 40), text_color="black")
        self.label_set.grid(row=0, column=0, padx=10, pady=10, sticky='w')

        self.sets = ["All", "Edges", "Nodes"]
        self.selected_set = ctk.StringVar()

        # Menu a tendina con il font personalizzato
        self.combobox_set = ctk.CTkComboBox(self.frame,
                                              values=self.sets,
                                              variable=self.selected_set,
                                              state='readonly',
                                              width=200,
                                              font=("helvetica", 30),
                                              dropdown_font=(
                                                  "helvetica", 30))  # Imposta il font per i valori nel menu a tendina
        self.combobox_set.grid(row=0, column=1, padx=10, pady=10)
        self.combobox_set.set(self.sets[0])

        self.process_button = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 40), width=200, height=40,
                                            fg_color="gray30")
        self.process_button.grid(row=1, column=0, padx=10, pady=20)

if __name__ == "__main__":
    root = ctk.CTk()
    app = CentralityAnalysis(root)
    root.mainloop()
