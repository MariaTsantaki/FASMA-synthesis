from unittest import TestCase
from FASMA.interpolation import solar_abundance


class TestSolarAbundance(TestCase):
    def test_return_values(self):
        element = 'Fe'
        index, abundance = solar_abundance(element)
        self.assertTrue(isinstance(index, int))
        self.assertTrue(isinstance(abundance, float))
        self.assertTrue(index == 26)
        self.assertTrue(abundance == 7.47)

    def test_none_return(self):
        element = 'tmp'
        index, abundance = solar_abundance(element)
        self.assertTrue(index is None)
        self.assertTrue(abundance is None)
