import tkinter as tk
from tkinter import messagebox
import requests
class GlobalSystemAccessPreview:
    def __init__(self):
        self.window = tk.Tk()
        self.window.title("Global System Access Preview")
        # Create widgets
        self.system_access_label = tk.Label(self.window, text="Global System Access Preview")
        self.system_access_label.pack()
        self.system_status = tk.StringVar()
        self.system_status.set("Offline")
        self.system_status_label = tk.Label(self.window, textvariable=self.system_status)
        self.system_status_label.pack()
        self.access_button = tk.Button(self.window, text="Check Access", command=self.check_access)
        self.access_button.pack()
    def check_access(self):
        try:
            response = requests.get("http://localhost:8080/access")
            if response.status_code == 200:
                self.system_status.set("Online")
            else:
                self.system_status.set("Offline")
        except requests.exceptions.RequestException as e:
            messagebox.showerror("Error", str(e))
    def run(self):
        self.window.mainloop()
if __name__ == "__main__":
    app = GlobalSystemAccessPreview()
    app.run()
```
