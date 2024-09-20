import customtkinter as ctk

class PathAnalysis:
    def __init__(self, root):
        self.root = root
        self.root.title("Path Analysis")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="DarkOliveGreen1")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per il posizionamento dei widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="DarkOliveGreen1")
        self.frame.pack(expand=True, padx=10, pady=10)

        self.label_start_list = ctk.CTkLabel(self.frame, text="Insert the start's list of residues:", font=("Helvetica", 40),
                                            text_color="black")
        self.label_start_list.grid(row=0, column=0, padx=10, pady=10)

        self.entry_start_list = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 30), justify="center")
        self.entry_start_list.grid(row=0, column=1, padx=10, pady=10)

        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30),
                                           command=self.confirm_start, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=0, column=2, padx=10, pady=10)

        self.label_end_list = ctk.CTkLabel(self.frame, text="Insert the end's list of residues:",
                                             font=("Helvetica", 40),
                                             text_color="black")
        self.label_end_list.grid(row=1, column=0, padx=10, pady=10)

        self.entry_end_list = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 30), justify="center")
        self.entry_end_list.grid(row=1, column=1, padx=10, pady=10)

        self.button_cutOff = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 30),
                                           command=self.confirm_end, fg_color="gray30", width=200)
        self.button_cutOff.grid(row=1, column=2, padx=10, pady=10)

        self.process_button = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), width=200, height=40,
                                            fg_color="gray30")
        self.process_button.grid(row=2, column=1, padx=10, pady=20)

    def confirm_start(self):
        self.entry_start_list.get()

    def confirm_end(self):
        self.entry_end_list.get()


if __name__ == "__main__":
    root = ctk.CTk()
    app = PathAnalysis(root)
    root.mainloop()
