/**
 * Test script to verify the backup flow works end-to-end
 * Run this after starting the backup server with: python backup_server.py
 */

const { backupToServer, buildBundle } = require('./dist/backup');

async function testBackupFlow() {
  console.log('🔵 Testing backup flow...\n');

  try {
    // Test 1: Build a bundle
    console.log('1. Building backup bundle...');
    const bundle = await buildBundle();
    console.log('   ✅ Bundle created:', {
      timestamp: bundle.timestamp,
      version: bundle.version,
      dataKeys: Object.keys(bundle.data)
    });

    // Test 2: Send backup to server
    console.log('\n2. Sending backup to server...');
    await backupToServer();
    console.log('   ✅ Backup sent successfully!');

    console.log('\n✅ All tests passed! Backup system is working correctly.');
    console.log('   Check the backups/ directory for saved backup files.');
    
  } catch (error) {
    console.error('\n❌ Error during backup:', error.message);
    process.exit(1);
  }
}

testBackupFlow();
