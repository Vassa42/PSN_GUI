import customtkinter as ctk
import customtkinter

from centrality_analysis import CentralityAnalysis
from connected_components import ConnectedComponents
from degree_analysis import DegreeAnalysis


class TopologyAnalysis:
    def __init__(self, root):
        self.root = root
        self.root.title("Topology Analysis")

        # Imposta il tema di customtkinter
        ctk.set_appearance_mode("dark")
        ctk.set_default_color_theme("blue")

        # Frame esterno con bordo
        self.outer_frame = ctk.CTkFrame(root, border_color="black", border_width=2, fg_color="PaleGreen1")
        self.outer_frame.pack(expand=True, fill='both', padx=10, pady=10)

        # Frame interno per il posizionamento dei widget
        self.frame = ctk.CTkFrame(self.outer_frame, fg_color="PaleGreen1")
        self.frame.pack(expand=True, padx=10, pady=10)

        self.label_metrics = customtkinter.CTkLabel(self.frame, text="Select the metrics:",
                                                   font=("Helvetica", 50), text_color="#000000")
        self.label_metrics.grid(row=0, column=0, padx=10, pady=10)

        # Variabili per i checkbox
        self.var_degree = customtkinter.IntVar()
        self.var_components = customtkinter.IntVar()
        self.var_centrality = customtkinter.IntVar()

        # Checkbox per i vari moduli
        self.check_degree = customtkinter.CTkCheckBox(self.frame, text="Degree Analysis", variable=self.var_degree,
                                                  font=("Helvetica", 40), text_color="#000000")
        self.check_components = customtkinter.CTkCheckBox(self.frame, text="Connected Components",
                                                  variable=self.var_components,
                                                  font=("Helvetica", 40), text_color="#000000")
        self.check_centrality = customtkinter.CTkCheckBox(self.frame, text="Centrality Analysis", variable=self.var_centrality,
                                                  font=("Helvetica", 40), text_color="#000000")

        self.check_degree.grid(row=1, column=0, sticky='w', padx=10)
        self.check_components.grid(row=2, column=0, sticky='w', padx=10)
        self.check_centrality.grid(row=3, column=0, sticky='w', padx=10)

        # Bottone la selezione
        self.button_select = customtkinter.CTkButton(self.frame, text="Select", command=self.open_metrics_analysis, font=("Helvetica", 40),
                                                     fg_color="#333333",  # Dark color for button
                                                     hover_color="#555555", width=200)  # Slightly lighter hover color
        self.button_select.grid(row=4, column=0, columnspan=3, pady=20)
    def open_metrics_analysis(self):
        if  self.var_degree.get() and self.var_components.get() and self.var_centrality.get():
            self.open_metrics(DegreeAnalysis)
            self.open_metrics(ConnectedComponents)
            self.open_metrics(CentralityAnalysis)

        elif self.var_degree.get() and self.var_components.get():
             self.open_metrics(DegreeAnalysis)
             self.open_metrics(ConnectedComponents)

        elif self.var_degree.get() and self.var_centrality.get():
             self.open_metrics(DegreeAnalysis)
             self.open_metrics(CentralityAnalysis)

        elif self.var_components.get() and self.var_centrality.get():
             self.open_metrics(ConnectedComponents)
             self.open_metrics(CentralityAnalysis)

        elif self.var_degree.get():
             self.open_metrics(DegreeAnalysis)

        elif self.var_components.get():
             self.open_metrics(ConnectedComponents)

        elif self.var_centrality.get():
             self.open_metrics(CentralityAnalysis)

    def open_metrics(self, metrics_class):
        new_window = customtkinter.CTkToplevel(self.root)
        metric = metrics_class(new_window)

if __name__ == "__main__":
    root = ctk.CTk()
    app = TopologyAnalysis(root)
    root.mainloop()
