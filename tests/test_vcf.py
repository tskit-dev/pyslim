"""
Test cases for VCF output in pyslim.
"""
import contextlib
import io
import tempfile

import numpy as np
import pytest
import vcf


@contextlib.contextmanager
def ts_to_pyvcf(ts, *args, **kwargs):
    """
    Returns a PyVCF reader for the specified tree sequence and arguments.
    """
    f = io.StringIO()
    ts.write_vcf(f, *args, **kwargs)
    f.seek(0)
    yield vcf.Reader(f)



class TestRoundTripIndividuals(PyslimTestCase):
    """
    Tests that we can round-trip genotype data through VCF using pyvcf.
    """

    def verify(self, ts):
        for indivs, _num_indivs in example_individuals(ts):
            with ts_to_pyvcf(ts, individuals=indivs) as vcf_reader:
                samples = []
                if indivs is None:
                    indivs = range(ts.num_individuals)
                for ind in map(ts.individual, indivs):
                    samples.extend(ind.nodes)
                for variant, vcf_row in itertools.zip_longest(
                    ts.variants(samples=samples), vcf_reader
                ):
                    assert vcf_row.POS == np.round(variant.site.position)
                    assert variant.alleles[0] == vcf_row.REF
                    assert list(variant.alleles[1:]) == vcf_row.ALT
                    j = 0
                    for individual, sample in itertools.zip_longest(
                        map(ts.individual, indivs), vcf_row.samples
                    ):
                        calls = sample.data.GT.split("|")
                        allele_calls = sample.gt_bases.split("|")
                        assert len(calls) == len(individual.nodes)
                        for allele_call, call in zip(allele_calls, calls):
                            assert int(call) == variant.genotypes[j]
                            assert allele_call == variant.alleles[variant.genotypes[j]]
                            j += 1

    def test_round_trip(self):
        for ts in self.get_slim_examples(nucleotides=True):
            self.verify(ts)

