#!/usr/bin/env bash
set -e

# THIS MUST BE RUN AS ROOT ON A LINUX SERVER WITH UFW INSTALLED.
# Adjust ports as needed for your setup.

echo "[*] Enabling UFW with default deny..."
ufw default deny incoming
ufw default allow outgoing

# Allow SSH (change 22 to your custom SSH port if you use one)
ufw allow 22/tcp

# Allow app port
APP_PORT=${APP_PORT:-8080}
ufw allow "${APP_PORT}/tcp"

# Log denied packets
ufw logging on

echo "[*] Enabling firewall..."
ufw --force enable

echo "[*] Firewall rules applied."
ufw status verbose
