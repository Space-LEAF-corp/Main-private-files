class DigitalRepellent:
    def __init__(self):
        self.cycles = 0
        self.threat_log = []
        self.firewall_rules = set()

    def scan_and_learn(self, traffic):
        for packet in traffic:
            if self.is_threat(packet):
                self.threat_log.append(packet)
                self.firewall_rules.add(self.generate_rule(packet))

    def is_threat(self, packet):
        return "malicious" in packet.get("tags", [])

    def generate_rule(self, packet):
        return f"BLOCK {packet['source_ip']}"

    def finalize_cycle(self):
        self.cycles += 1
        if self.cycles % 3 == 0:
            print("🔒 Finalizing hardening phase...")
            self.deploy_rules()

    def deploy_rules(self):
        print(f"Deploying {len(self.firewall_rules)} rules.")
        # Simulate rule deployment
        for rule in self.firewall_rules:
            print(rule)

# Simulated usage
repellent = DigitalRepellent()
for day in range(1, 366):
    traffic = [{"source_ip": f"192.168.0.{day%255}", "tags": ["malicious"] if day % 5 == 0 else []}]
    repellent.scan_and_learn(traffic)
    repellent.finalize_cycle()
