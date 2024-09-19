import customtkinter as ctk
from customtkinter import filedialog
import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt

class PreliminaryPlot:
    def __init__(self, root):
        self.root = root
        self.root.title("Preliminary Plot")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="light goldenrod")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per il posizionamento dei widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="light goldenrod")
        self.frame.pack(expand=True, padx=10, pady=10)

        # Valori Min, Max e Step per il filtraggio
        self.label_min = ctk.CTkLabel(self.frame, text="Minimum value (%):", font=("Helvetica", 50), text_color="black")
        self.label_min.grid(row=0, column=0, padx=10, pady=10, sticky='w')

        self.entry_min = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 40))
        self.entry_min.grid(row=0, column=1, padx=10, pady=10)
        self.entry_min.insert(0, "0.0")

        self.button_min = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40), command=self.update_min,
                                        width=200, height=40, fg_color="gray30")
        self.button_min.grid(row=0, column=2, padx=10, pady=10)

        self.label_max = ctk.CTkLabel(self.frame, text="Maximum value (%):", font=("Helvetica", 50), text_color="black")
        self.label_max.grid(row=1, column=0, padx=10, pady=10, sticky='w')

        self.entry_max = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 40))
        self.entry_max.grid(row=1, column=1, padx=10, pady=10)
        self.entry_max.insert(0, "100.0")

        self.button_max = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40), command=self.update_max,
                                        width=200, height=40, fg_color="gray30")
        self.button_max.grid(row=1, column=2, padx=10, pady=10)

        self.label_step = ctk.CTkLabel(self.frame, text="Step value (%):", font=("Helvetica", 50), text_color="black")
        self.label_step.grid(row=2, column=0, padx=10, pady=10, sticky='w')

        self.entry_step = ctk.CTkEntry(self.frame, width=200, height=40, font=("Helvetica", 40))
        self.entry_step.grid(row=2, column=1, padx=10, pady=10)
        self.entry_step.insert(0, "1.0")

        self.button_step = ctk.CTkButton(self.frame, text="Confirm", font=("Helvetica", 40), command=self.update_step,
                                         width=200, height=40, fg_color="gray30")
        self.button_step.grid(row=2, column=2, padx=10, pady=10)

        self.process_button = ctk.CTkButton(self.frame, text="Process", font=("Helvetica", 40), width=200, height=40, fg_color="gray30")
        self.process_button.grid(row=3, column=1, padx=10, pady=20)
    def update_min(self):
        try:
            self.min_value = float(self.entry_min.get())
        except ValueError:
            print("Error: Invalid min value.")

    def update_max(self):
        try:
            self.max_value = float(self.entry_max.get())
        except ValueError:
            print("Error: Invalid max value.")

    def update_step(self):
        try:
            self.step_value = float(self.entry_step.get())
        except ValueError:
            print("Error: Invalid step value.")


if __name__ == "__main__":
    root = ctk.CTk()
    app = PreliminaryPlot(root)
    root.mainloop()
