class ThermalHeatGenerator:
    def __init__(self):
        self.last_power_output_w = 0.0
        self.last_heat_flux = 0.0

    def sense_heat(self, surface_temp, inner_temp):
        # Simple gradient proxy
        gradient = max(0.0, surface_temp - inner_temp)
        self.last_heat_flux = gradient
        return gradient

    def generate_power(self, gradient):
        # Placeholder: more gradient = more emergency power
        efficiency = 0.01  # conceptual
        self.last_power_output_w = gradient * efficiency
        return self.last_power_output_w
