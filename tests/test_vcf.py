"""
Test cases for VCF output in pyslim.
"""

import pyslim
import contextlib
import io
import tempfile
import numpy as np
import itertools
import tests
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


def example_individuals(ts, ploidy=1):
    if ts.num_individuals == 0:
        yield None, ts.num_samples / ploidy
    else:
        yield None, ts.num_individuals
        yield list(range(ts.num_individuals)), ts.num_individuals
    if ts.num_individuals > 3:
        n = ts.num_individuals - 2
        yield list(range(n)), n
        yield 2 + np.random.choice(np.arange(n), n, replace=False), n


class TestNucleotideBased(tests.PyslimTestCase):
    """
    Tests that we can round-trip genotype data through VCF using pyvcf
    for nucleotide models.
    """

    def verify(self, ts):
        for indivs, _num_indivs in example_individuals(ts):
            with ts_to_pyvcf(ts, individuals=indivs, nucleotide_based=True) as vcf_reader:
                samples = []
                if indivs is None:
                    indivs = range(ts.num_individuals)
                for ind in map(ts.individual, indivs):
                    samples.extend(ind.nodes)
                for variant, vcf_row in itertools.zip_longest(
                    ts.variants(samples=samples), vcf_reader
                ):
                    assert vcf_row.POS == np.round(variant.site.position)
                    assert ts.reference_sequence[int(variant.site.position)] == vcf_row.REF
                    alleles = [ts.reference_sequence[int(variant.site.position)]] 
                    alleles += ["".join(
                            [pyslim.NUCLEOTIDES[u['nucleotide']] for u in m.metadata['mutation_list']]
                        ) for m in variant.site.mutations]
                    vcf_alleles = [vcf_row.REF] + vcf_row.ALT
                    j = 0
                    for individual, sample in itertools.zip_longest(
                        map(ts.individual, indivs), vcf_row.samples
                    ):
                        calls = sample.data.GT.split("|")
                        allele_calls = sample.gt_bases.split("|")
                        assert len(calls) == len(individual.nodes)
                        for allele_call, call in zip(allele_calls, calls):
                            assert vcf_alleles[int(call)] == alleles[variant.genotypes[j]]
                            assert allele_call == alleles[variant.genotypes[j]]
                            j += 1

    def test_round_trip(self):
        for ts in self.get_slim_examples(nucleotides=True):
            self.verify(ts)


class TestNonNucleotideBased(tests.PyslimTestCase):
    """
    Tests that we can round-trip genotype data through VCF using pyvcf
    for non-nucleotide models.
    """

    def verify(self, ts):
        for indivs, _num_indivs in example_individuals(ts):
            with ts_to_pyvcf(ts, individuals=indivs, nucleotide_based=False) as vcf_reader:
                samples = []
                if indivs is None:
                    indivs = range(ts.num_individuals)
                for ind in map(ts.individual, indivs):
                    samples.extend(ind.nodes)
                for variant, vcf_row in itertools.zip_longest(
                    ts.variants(samples=samples), vcf_reader
                ):
                    assert vcf_row.POS == np.round(variant.site.position)
                    assert pyslim.NUCLEOTIDES[0] == vcf_row.REF
                    vcf_alleles = [vcf_row.REF] + vcf_row.ALT
                    vcf_genotypes = []
                    for a in vcf_alleles:
                        n = 0
                        b = 1
                        for x in str(a):
                            n += pyslim.NUCLEOTIDES.index(x) * b
                            b *= 4
                        vcf_genotypes.append(n)
                    j = 0
                    for individual, sample in itertools.zip_longest(
                        map(ts.individual, indivs), vcf_row.samples
                    ):
                        calls = sample.data.GT.split("|")
                        allele_calls = sample.gt_bases.split("|")
                        assert len(calls) == len(individual.nodes)
                        for allele_call, call in zip(allele_calls, calls):
                            assert vcf_genotypes[int(call)] == variant.genotypes[j]
                            j += 1

    def test_round_trip(self):
        ts = self.get_slim_example(name="recipe_WF")
        self.verify(ts)

