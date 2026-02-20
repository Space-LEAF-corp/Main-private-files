"""
CLI Interface
-------------
Command-line interface for submitting GQI&EA inquiries.
"""

from src.core.tag_generator import generate_user_tag
from src.core.timestamp_engine import get_timestamp
from src.core.routing_engine import route_inquiry

def submit_from_cli(username: str, question: str):
    tag = generate_user_tag(username)
    timestamp = get_timestamp()
    route = route_inquiry(question)

    entry = (
        f"GQI&EA | user: {username} | {timestamp} | tag: {tag}\n"
        f"{question}\n"
        f"Routed to: {route}"
    )

    print(entry)
