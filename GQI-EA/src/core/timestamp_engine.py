from datetime import datetime

def get_timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M %Z")
