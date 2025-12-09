# Backup System

This directory contains the backup functionality for the Main Private Files system, providing both client-side and server-side backup capabilities.

## Overview

The backup system consists of:
- **TypeScript/JavaScript Client**: Functions to build and send backup bundles
- **Python Backup Server**: HTTP server to receive and store backups on port 8787

## Client-Side (TypeScript)

### Installation

```bash
npm install
npm run build
```

### Usage

```typescript
import { backupToServer, buildBundle } from './src/backup';

// Build a backup bundle
const bundle = await buildBundle();
console.log(bundle);

// Send backup to server
await backupToServer();
```

### API

#### `buildBundle(): Promise<BackupBundle>`

Creates a bundle of data to be backed up. The bundle includes:
- `timestamp`: ISO timestamp of when the bundle was created
- `version`: Version of the backup format
- `data`: Object containing users, sessions, firewall state, and metadata

#### `backupToServer(): Promise<void>`

Sends the backup bundle to the backup server at `http://localhost:8787/backup` via POST request.

## Server-Side (Python)

### Running the Server

```bash
python backup_server.py
```

The server will start on `http://localhost:8787` and listen for POST requests to `/backup`.

### Backup Endpoint

**POST /backup**

Accepts a JSON backup bundle and saves it to the `backups/` directory.

**Request Body:**
```json
{
  "timestamp": "2024-01-01T00:00:00.000Z",
  "version": "1.0.0",
  "data": {
    "users": [],
    "sessions": [],
    "firewall_state": {},
    "metadata": {}
  }
}
```

**Response:**
```json
{
  "status": "success",
  "message": "Backup saved successfully",
  "filename": "backup_20240101_120000.json",
  "timestamp": "20240101_120000",
  "size": 256
}
```

### Backup Storage

Backups are stored in the `backups/` directory with filenames in the format:
```
backup_YYYYMMDD_HHMMSS.json
```

## Testing

### TypeScript Tests

```bash
npm run build
npm test
```

### Python Tests

```bash
python -m pytest tests/test_backup_server.py
# or
python -m unittest tests.test_backup_server
```

## Integration Example

1. Start the backup server:
```bash
python backup_server.py
```

2. In another terminal, build and run the TypeScript client:
```bash
npm install
npm run build
node -e "require('./dist/backup').backupToServer().then(() => console.log('Backup completed'))"
```

3. Check the `backups/` directory for the saved backup file.

## Security Considerations

- The backup server currently runs on localhost by default for security
- In production, consider adding authentication to the backup endpoint
- Backup files should be encrypted if they contain sensitive data
- Use HTTPS instead of HTTP for production deployments

## Future Enhancements

- Add encryption for backup data
- Implement backup authentication/authorization
- Add backup scheduling capabilities
- Implement backup rotation and cleanup
- Add backup verification and integrity checks
