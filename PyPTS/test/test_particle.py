import unittest
from PyPTS import particle

class TestParticleData(unittest.TestCase):

    def test_init(self):
        pdata = particle.ParticleData(10)

        # re-initialized
        pdata = particle.ParticleData(100)

    def test_backup(self):
        pdata = particle.ParticleData("backup_100.pdata")
        pdata = particle.ParticleData("backup_100.pdata",copy=True)
        self.assertEqual(pdata.n_, 10)


if __name__ == "__main__":
    unittest.main()