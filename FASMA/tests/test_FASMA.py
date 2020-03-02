from unittest import TestCase
from FASMA.interpolation import solar_abundance
from FASMA import FASMA


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


class TestAbundances(TestCase):
    def test_minimization_values(self):
        elements = ['Li', 'Na', 'Mg', 'Al', 'Si', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Ni']
        for el in elements:
            options = {'linelist':'../rawLinelist/elements.lst',
            'intervals_file':'../rawLinelist/intervals_elements.lst', 'element':el,
            'resolution':115000, 'observations': '../spectra/Sun_HARPS.fits'}
            r = FASMA(**options)
            result = r.result()
            self.assertTrue(result['abund'] is not None)

    def test_minimization(self):
        options = {'linelist':'../rawLinelist/linelist.lst',
        'intervals_file':'../rawLinelist/intervals.lst',
        'resolution':115000, 'observations': '../spectra/Sun_HARPS.fits'}
        r = FASMA(**options)
        result = r.result()
        self.assertTrue(result is None)


if __name__ == '__main__':
    p = TestSolarAbundance()
    p.test_return_values()
    p.test_none_return()
    r = TestAbundances()
    r.test_minimization_values()
    r.test_minimization()
