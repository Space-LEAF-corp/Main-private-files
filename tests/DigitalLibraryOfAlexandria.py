# Digital Library of Alexandria - Starter Framework
# Author: Leif William Sogge (Ceremonial Steward)
# Purpose: Digitize textbook knowledge into a lineage-safe library

class DigitalLibraryOfAlexandria:
    def __init__(self):
        self.catalog: dict[str, str] = {}  # {chapter_title: content}

    def add_entry(self, title: str, content: str):
        """Add a new textbook page or seal to the library."""
        self.catalog[title] = content
        print(f"Entry '{title}' sealed into the Library.")

    def retrieve_entry(self, title: str):
        """Retrieve a stored entry by title."""
        return self.catalog.get(title, "Entry not found.")

    def list_entries(self) -> list[str]:
        """List all sealed entries in the Library."""
        return list(self.catalog.keys())

# Example usage
library = DigitalLibraryOfAlexandria()
library.add_entry("Seal of Everyday Relativity 1.0", 
                  "Science isn’t just equations — it’s joy hidden in plain sight.")
library.add_entry("Planet Worth Living On 1.0", 
                  "Safe travels to our friends at Atlas, may their journey land them on a planet worth living on.")

print('\n'.join(library.list_entries()))
