"""CLI to run the dual-layer server."""
import argparse
import sys
import time

from server import DualLayerServer


def main(argv=None):
    parser = argparse.ArgumentParser(prog="firewall-server")
    parser.add_argument(
        "--host",
        default="localhost",
        help="Bind address (default: localhost)",
    )
    parser.add_argument(
        "--socket-port",
        type=int,
        default=9000,
        help="Socket server port (default: 9000)",
    )
    parser.add_argument(
        "--http-port",
        type=int,
        default=9001,
        help="HTTP server port (default: 9001)",
    )

    args = parser.parse_args(argv)

    server = DualLayerServer(
        host=args.host,
        socket_port=args.socket_port,
        http_port=args.http_port,
    )
    server.start()

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\nShutting down...")
        server.stop()


if __name__ == "__main__":
    main(sys.argv[1:])
