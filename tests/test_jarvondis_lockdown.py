import os
import unittest
from main import load_policy
from jarvondis_lockdown import LockdownPolicy, AdministrativeLockdown

class TestLockdownPolicyEnvSecret(unittest.TestCase):
    def setUp(self):
        # Ensure env var unset for tests that expect failure
        self.orig = os.environ.pop('JARVONDIS_ADMIN_SECRET', None)

    def tearDown(self):
        # Restore env
        if self.orig is not None:
            os.environ['JARVONDIS_ADMIN_SECRET'] = self.orig
        else:
            os.environ.pop('JARVONDIS_ADMIN_SECRET', None)

    def test_load_policy_requires_env(self):
        # jarvondis_policy.json on branch requires admin_secret via env
        with self.assertRaises(ValueError):
            load_policy('jarvondis_policy.json')

    def test_load_policy_with_env_and_respond(self):
        os.environ['JARVONDIS_ADMIN_SECRET'] = 'test_env_secret_123'
        policy = load_policy('jarvondis_policy.json')
        self.assertEqual(policy.admin_secret, 'test_env_secret_123')
        guard = AdministrativeLockdown(policy)

        # Allowlisted client
        r = guard.respond('dimitri', 'Dimitri/1.0', 'status')
        self.assertEqual(r['status'], 'ok')
        self.assertEqual(r['mode'], 'minimal_lockdown')
        self.assertEqual(r['owner'], policy.owner_id)

        # Non-allowlisted client
        r2 = guard.respond('miko', 'Miko/1.0', 'status')
        self.assertEqual(r2['status'], 'blocked')

if __name__ == '__main__':
    unittest.main()
