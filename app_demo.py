from fastapi import FastAPI, Request, HTTPException
from pydantic import BaseModel
from datetime import datetime
import uuid

# ============================================================
#  BASIC UTILS (HASHING, ENCRYPTION, CODES)
# ============================================================


def hash_value(value: str) -> str:
	return f"hash({value})"


def verify_hash(value: str, hashed: str) -> bool:
	return hashed == f"hash({value})"


def generate_code(length=6) -> str:
	return "123456"  # Replace with real generator


def encrypt_blob(data: str) -> str:
	return f"enc({data})"


def decrypt_blob(data: str) -> str:
	if data.startswith("enc(") and data.endswith(")"):
		return data[4:-1]
	return data


def uuid_str() -> str:
	return str(uuid.uuid4())


# ============================================================
#  MOCK EXTERNAL SERVICES
# ============================================================


async def send_email(email: str, subject: str, body: str):
	print(f"[EMAIL] To:{email} | {subject} | {body}")


async def lookup_geo(ip: str) -> str:
	return "US-FL"


async def verify_satellite_uplink(token: str, context=None):
	return token == "VALID_UPLINK"


# ============================================================
#  IN-MEMORY DATABASES
# ============================================================

DB_USERS = {}
DB_EMAIL_CODES = {}
DB_IDENTITY = {}
DB_GUARDIANS = {}
DB_GLOBES = {}
DB_EVENTS = {}
DB_SESSIONS = {}

# ============================================================
#  REQUEST MODELS
# ============================================================


class SignupStart(BaseModel):
	email: str
	username: str
	password: str


class VerifyEmail(BaseModel):
	userId: str
	code: str


class VerifyGovID(BaseModel):
	userId: str
	govIdData: str
	satelliteUplinkToken: str | None = None


class SetupGuardian(BaseModel):
	userId: str
	spiritAnimal: str
	loginImageUrl: str


class LoginStart(BaseModel):
	email: str
	password: str


class LoginComplete(BaseModel):
	userId: str
	code: str


# ============================================================
#  FASTAPI APP
# ============================================================

app = FastAPI()


# ============================================================
#  DB HELPERS
# ============================================================
async def create_user(data):
	user_id = uuid_str()
	DB_USERS[user_id] = {
	 "email": data["email"],
	 "username": data["username"],
	 "passwordHash": data["passwordHash"],
	 "emailVerified": False
	}
	return user_id


async def save_email_code(user_id, code):
	DB_EMAIL_CODES[user_id] = code


async def verify_email_code(user_id, code):
	return DB_EMAIL_CODES.get(user_id) == code


async def mark_email_verified(user_id):
	DB_USERS[user_id]["emailVerified"] = True


async def get_user_email(user_id):
	return DB_USERS[user_id]["email"]


async def save_identity_token(token):
	DB_IDENTITY[token["userId"]] = token


async def get_identity_token(user_id):
	return DB_IDENTITY.get(user_id)


async def save_guardian_profile(profile):
	DB_GUARDIANS[profile["userId"]] = profile


async def get_guardian_profile(user_id):
	return DB_GUARDIANS.get(user_id)


async def create_empty_gyro_globe(globe):
	DB_GLOBES[globe["userId"]] = globe


async def get_gyro_globe(user_id):
	return DB_GLOBES.get(user_id)


async def log_gyro_event(event):
	DB_EVENTS.setdefault(event["userId"], []).append(event)


async def update_gyro_globe(user_id):
	events = DB_EVENTS.get(user_id, [])
	DB_GLOBES[user_id]["visualizationState"] = f"nodes:{len(events)}"
	DB_GLOBES[user_id]["lastUpdatedAt"] = datetime.utcnow()


async def evaluate_risk(session, guardian):
	return "LOW"  # Placeholder


async def sponge_pass(user_id):
	return True  # Placeholder


async def create_session_token(session):
	token = uuid_str()
	DB_SESSIONS[token] = session
	return token


async def destroy_session(token):
	DB_SESSIONS.pop(token, None)


# ============================================================
#  ROUTES
# ============================================================


@app.post("/signup/start")
async def signup_start(data: SignupStart):
	user_id = await create_user(
	 {
	  "email": data.email,
	  "username": data.username,
	  "passwordHash": hash_value(data.password)
	 }
	)

	code = generate_code()
	await save_email_code(user_id, code)
	await send_email(data.email, "Your Verification Key", f"Key: {code}")

	return {"userId": user_id, "status": "EMAIL_KEY_SENT"}


@app.post("/signup/verify-email")
async def signup_verify_email(data: VerifyEmail):
	if not await verify_email_code(data.userId, data.code):
		raise HTTPException(400, "INVALID_KEY_OR_EXPIRED")

	await mark_email_verified(data.userId)
	return {"status": "EMAIL_VERIFIED"}


@app.post("/signup/verify-gov-id")
async def signup_verify_gov_id(data: VerifyGovID):
	gov_hash = hash_value(data.govIdData)

	token = {
	 "userId": data.userId,
	 "email": await get_user_email(data.userId),
	 "govIdHash": gov_hash,
	 "createdAt": datetime.utcnow(),
	 "lastVerifiedAt": datetime.utcnow()
	}

	if data.satelliteUplinkToken:
		verified = await verify_satellite_uplink(data.satelliteUplinkToken)
		if not verified:
			raise HTTPException(400, "INVALID_SATELLITE_UPLINK")
		token["satelliteUplink"] = encrypt_blob(data.satelliteUplinkToken)

	await save_identity_token(token)
	return {
	 "status": "IDENTITY_VERIFIED",
	 "satelliteLinked": "satelliteUplink" in token
	}


@app.post("/signup/setup-guardian")
async def signup_setup_guardian(data: SetupGuardian):
	guardian = {
	 "userId": data.userId,
	 "spiritAnimal": data.spiritAnimal,
	 "loginImageUrl": data.loginImageUrl,
	 "securityPreset": "AGGRESSIVE",
	 "firewallActive": True
	}

	await save_guardian_profile(guardian)

	globe = {
	 "userId": data.userId,
	 "globeId": uuid_str(),
	 "encryptedEventsBlob": encrypt_blob("[]"),
	 "visualizationState": f"theme:{data.spiritAnimal}",
	 "lastUpdatedAt": datetime.utcnow()
	}

	await create_empty_gyro_globe(globe)

	return {"status": "GUARDIAN_CONFIGURED", "firewall": "ACTIVE"}


@app.post("/login/start")
async def login_start(data: LoginStart):
	user = next((u for u in DB_USERS.items() if u[1]["email"] == data.email), None)
	if not user:
		raise HTTPException(401, "INVALID_CREDENTIALS")

	user_id, user_data = user

	if not verify_hash(data.password, user_data["passwordHash"]):
		raise HTTPException(401, "INVALID_CREDENTIALS")

	identity = await get_identity_token(user_id)
	if not identity:
		raise HTTPException(403, "UNVERIFIED_IDENTITY")

	code = generate_code()
	await save_email_code(user_id, code)
	await send_email(data.email, "Your Login Key", f"Key: {code}")

	return {"userId": user_id, "status": "EMAIL_KEY_SENT"}


@app.post("/login/complete")
async def login_complete(data: LoginComplete, request: Request):
	if not await verify_email_code(data.userId, data.code):
		raise HTTPException(400, "INVALID_KEY")

	session = {
	 "userId": data.userId,
	 "ip": request.client.host,
	 "userAgent": request.headers.get("user-agent", ""),
	 "geo": await lookup_geo(request.client.host)
	}

	guardian = await get_guardian_profile(data.userId)
	if not guardian or not guardian["firewallActive"]:
		raise HTTPException(403, "GUARDIAN_FIREWALL_INACTIVE")

	risk = await evaluate_risk(session, guardian)
	if risk == "HIGH":
		raise HTTPException(403, "GUARDIAN_FIREWALL_BLOCK")

	identity = await get_identity_token(data.userId)
	if "satelliteUplink" in identity:
		uplink = decrypt_blob(identity["satelliteUplink"])
		ok = await verify_satellite_uplink(uplink, {"ip": session["ip"]})
		if not ok:
			raise HTTPException(403, "SATELLITE_UPLINK_FAILURE")

		await log_gyro_event(
		 {
		  "userId": data.userId,
		  "eventId": uuid_str(),
		  "type": "SATELLITE_UPLINK",
		  "timestamp": datetime.utcnow(),
		  "metadata": {
		   "status": "VERIFIED"
		  }
		 }
		)

	await log_gyro_event(
	 {
	  "userId": data.userId,
	  "eventId": uuid_str(),
	  "type": "LOGIN",
	  "timestamp": datetime.utcnow(),
	  "metadata": session
	 }
	)

	await update_gyro_globe(data.userId)
	await sponge_pass(data.userId)

	token = await create_session_token(session)

	return {
	 "token": token,
	 "guardian": guardian,
	 "status": "LOGIN_OK",
	 "planetGlobe": await get_gyro_globe(data.userId)
	}

