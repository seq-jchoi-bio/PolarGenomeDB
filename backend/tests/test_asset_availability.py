import tempfile
from pathlib import Path
import unittest

from backend.routes.assets import asset_record, prefix_asset_record


class AssetAvailabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.prefix = Path(self.temp.name) / 'genome'

    def files(self, suffixes):
        for suffix in suffixes:
            Path(str(self.prefix) + suffix).write_bytes(b'nonempty-test-fixture')

    def ready(self, kind):
        return prefix_asset_record(str(self.prefix), {'asset_type': kind})['exists']

    def test_empty_and_directory_are_not_files(self):
        self.prefix.touch()
        self.assertFalse(asset_record(str(self.prefix))['exists'])
        self.assertFalse(asset_record(self.temp.name)['exists'])

    def test_bowtie_requires_all_six_components(self):
        self.files(['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2'])
        self.assertFalse(self.ready('bowtie2_genome'))
        self.files(['.rev.2.bt2'])
        self.assertTrue(self.ready('bowtie2_genome'))

    def test_large_bowtie(self):
        self.files([f'.{part}.bt2l' for part in ('1', '2', '3', '4', 'rev.1', 'rev.2')])
        self.assertTrue(self.ready('bowtie2_genome'))

    def test_blast_core_files(self):
        self.files(['.nhr', '.nin'])
        self.assertFalse(self.ready('blast_genome'))
        self.files(['.nsq'])
        self.assertTrue(self.ready('blast_genome'))

    def test_split_blast(self):
        self.files(['.00.phr', '.00.pin', '.00.psq', '.01.psq'])
        self.assertFalse(self.ready('blast_protein'))
        self.files(['.01.phr', '.01.pin'])
        self.assertTrue(self.ready('blast_protein'))

    def test_diamond_log_does_not_count(self):
        self.files(['.dmnd.log'])
        self.assertFalse(self.ready('diamond_protein'))
        self.files(['.dmnd'])
        self.assertTrue(self.ready('diamond_protein'))


if __name__ == '__main__':
    unittest.main()
