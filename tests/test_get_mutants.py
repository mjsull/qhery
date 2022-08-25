from qhery.get_mutants import mutantFinder
import pytest

working_dir='data'
sample_name='example'

@pytest.fixture()
def mf():
    return mutantFinder(working_dir=working_dir, sample_name='example')

def test_mutantFinder_init():
    mf = mutantFinder(working_dir=working_dir, sample_name='example')
    assert mf.bcftools_vcf == f'{working_dir}/{sample_name}.bcftools.vcf'
    assert mf.csq_file == f'{working_dir}/{sample_name}.csq.tsv'

def test_parse_csq(mf: mutantFinder):
    res = mf.parse_csq()
    res = []