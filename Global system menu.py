# Import necessary libraries
import tkinter as tk
class GlobalSystemMenu:
    def __init__(self):
        self.window = tk.Tk()
        self.window.title("Global System Menu")
        # Create widgets
        self.system_status_label = tk.Label(self.window, text="System Status:")
        self.system_status_label.pack()
        self.system_status_value = tk.Label(self.window, text="Offline")
        self.system_status_value.pack()
        self.system_info_label = tk.Label(self.window, text="System Information:")
        self.system_info_label.pack()
        self.system_info_value = tk.Label(self.window, text="Not Available")
        self.system_info_value.pack()
        self.system_logs_label = tk.Label(self.window, text="System Logs:")
        self.system_logs_label.pack()
        self.system_logs_value = tk.Label(self.window, text="Not Available")
        self.system_logs_value.pack()
        self.system_settings_label = tk.Label(self.window, text="System Settings:")
        self.system_settings_label.pack()
        self.system_settings_value = tk.Label(self.window, text="Not Available")
        self.system_settings_value.pack()
    def run(self):
        self.window.mainloop()
if __name__ == "__main__":
    menu = GlobalSystemMenu()
    menu.run()
```