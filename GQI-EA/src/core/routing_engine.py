"""
Routing Engine
--------------
Determines where an inquiry should be sent based on its content.
"""

def route_inquiry(inquiry: str) -> str:
    """
    Placeholder routing logic.
    """
    if "reset" in inquiry.lower():
        return "technical-support"
    if "where" in inquiry.lower():
        return "general-info"
    return "general-support"
