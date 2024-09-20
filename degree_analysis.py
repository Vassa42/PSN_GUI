import customtkinter as ctk


class DegreeAnalysis:
    def __init__(self, root):
        self.root = root
        self.root.title("Degree Analysis")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="DarkSeaGreen1")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per il posizionamento dei widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="DarkSeaGreen1")
        self.frame.pack(expand=True, padx=10, pady=10)

        self.process_button = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 30), width=200, height=40,
                                            fg_color="gray30")
        self.process_button.grid(row=0, column=0, padx=10, pady=20)

if __name__ == "__main__":
    root = ctk.CTk()
    app = DegreeAnalysis(root)
    root.mainloop()
