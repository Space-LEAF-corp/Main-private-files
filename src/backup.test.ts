/**
 * Tests for backup module
 */

import { buildBundle, backupToServer } from './backup';
import * as assert from 'assert';
import { test, describe } from 'node:test';

describe('Backup Module', () => {
  test('buildBundle creates a valid bundle', async () => {
    const bundle = await buildBundle();
    
    // Verify bundle structure
    assert.ok(bundle, 'Bundle should exist');
    assert.ok(bundle.timestamp, 'Bundle should have a timestamp');
    assert.ok(bundle.version, 'Bundle should have a version');
    assert.ok(bundle.data, 'Bundle should have data');
    
    // Verify data structure
    assert.ok(Array.isArray(bundle.data.users), 'Users should be an array');
    assert.ok(Array.isArray(bundle.data.sessions), 'Sessions should be an array');
    assert.ok(typeof bundle.data.firewall_state === 'object', 'Firewall state should be an object');
    assert.ok(typeof bundle.data.metadata === 'object', 'Metadata should be an object');
    
    console.log('✓ buildBundle test passed');
  });

  test('buildBundle includes metadata', async () => {
    const bundle = await buildBundle();
    
    assert.ok(bundle.data.metadata, 'Metadata should exist');
    assert.strictEqual(bundle.data.metadata.source, 'main-private-files', 'Source should be set');
    assert.ok(typeof bundle.data.metadata.buildTime === 'number', 'Build time should be a number');
    
    console.log('✓ buildBundle metadata test passed');
  });

  test('buildBundle timestamp is valid ISO string', async () => {
    const bundle = await buildBundle();
    const timestamp = new Date(bundle.timestamp);
    
    assert.ok(!isNaN(timestamp.getTime()), 'Timestamp should be valid');
    assert.ok(timestamp.getTime() > 0, 'Timestamp should be positive');
    
    console.log('✓ buildBundle timestamp test passed');
  });
});
