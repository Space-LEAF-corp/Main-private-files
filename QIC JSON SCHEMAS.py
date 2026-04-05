import json
import time
from dataclasses import dataclass, asdict
from jsonschema import validate

# ---------------------------------------------------------
# QIC JSON Schemas
# ---------------------------------------------------------

SCHEMA_QIC_DNA = {
    "type": "object",
    "required": ["type", "version", "identity_hash", "public_key", "timestamp"],
    "properties": {
        "type": {"const": "QIC-DNA"},
        "version": {"type": "string"},
        "identity_hash": {"type": "string"},
        "public_key": {"type": "string"},
        "timestamp": {"type": "integer"},
        "metadata": {"type": "object"}
    }
}

SCHEMA_QIC_SYS = {
    "type": "object",
    "required": ["type", "version", "system_id", "recovery_key", "nonce", "timestamp"],
    "properties": {
        "type": {"const": "QIC-SYS"},
        "version": {"type": "string"},
        "system_id": {"type": "string"},
        "recovery_key": {"type": "string"},
        "nonce": {"type": "string"},
        "timestamp": {"type": "integer"},
        "signature": {"type": "string"}
    }
}

SCHEMA_QIC_USER = {
    "type": "object",
    "required": ["type", "version", "user_id", "recovery_token", "timestamp"],
    "properties": {
        "type": {"const": "QIC-USER"},
        "version": {"type": "string"},
        "user_id": {"type": "string"},
        "recovery_token": {"type": "string"},
        "timestamp": {"type": "integer"},
        "cooldown": {"type": "integer"},
        "signature": {"type": "string"}
    }
}

SCHEMA_QIC_AUTH = {
    "type": "object",
    "required": ["type", "version", "credential_id", "public_key", "timestamp"],
    "properties": {
        "type": {"const": "QIC-AUTH"},
        "version": {"type": "string"},
        "credential_id": {"type": "string"},
        "public_key": {"type": "string"},
        "timestamp": {"type": "integer"},
        "binding": {"type": "string"},
        "signature": {"type": "string"}
    }
}

SCHEMA_QIC_OVR = {
    "type": "object",
    "required": [
        "type", "version", "override_id", "scope",
        "one_time_token", "timestamp", "expires"
    ],
    "properties": {
        "type": {"const": "QIC-OVR"},
        "version": {"type": "string"},
        "override_id": {"type": "string"},
        "scope": {"type": "string"},
        "one_time_token": {"type": "string"},
        "timestamp": {"type": "integer"},
        "expires": {"type": "integer"},
        "dual_control_required": {"type": "boolean"},
        "signature": {"type": "string"}
    }
}

# ---------------------------------------------------------
# Dataclasses for each QIC payload
# ---------------------------------------------------------

@dataclass
class QICBase:
    type: str
    version: str = "1.0"
    timestamp: int = int(time.time())

    def to_json(self, schema):
        payload = asdict(self)
        validate(payload, schema)
        return json.dumps(payload, separators=(",", ":"))


@dataclass
class QIC_DNA(QICBase):
    identity_hash: str = ""
    public_key: str = ""
    metadata: dict = None

    def to_json(self):
        return super().to_json(SCHEMA_QIC_DNA)


@dataclass
class QIC_SYS(QICBase):
    system_id: str = ""
    recovery_key: str = ""
    nonce: str = ""
    signature: str = ""

    def to_json(self):
        return super().to_json(SCHEMA_QIC_SYS)


@dataclass
class QIC_USER(QICBase):
    user_id: str = ""
    recovery_token: str = ""
    cooldown: int = 0
    signature: str = ""

    def to_json(self):
        return super().to_json(SCHEMA_QIC_USER)


@dataclass
class QIC_AUTH(QICBase):
    credential_id: str = ""
    public_key: str = ""
    binding: str = ""
    signature: str = ""

    def to_json(self):
        return super().to_json(SCHEMA_QIC_AUTH)


@dataclass
class QIC_OVR(QICBase):
    override_id: str = ""
    scope: str = ""
    one_time_token: str = ""
    expires: int = 0
    dual_control_required: bool = True
    signature: str = ""

    def to_json(self):
        return super().to_json(SCHEMA_QIC_OVR)