from typing import TypedDict

class ImmutableLabel(TypedDict):
    name: str
    version: str
    purpose: str
    harm_prohibition: bool
    repair_only: bool
    warnings_target: str
    covenant: str

immutable_label: ImmutableLabel = {
    "name": "Federated Stewardship Runtime",
    "version": "1.4",
    "purpose": "Simultaneous improvement across teams via repair-only governance and transparent labeling",
    "harm_prohibition": True,
    "repair_only": True,
    "warnings_target": "intruders",
    # Optional: covenant reference
    "covenant": "Seal of Immutable Labeling v1.0"
}
