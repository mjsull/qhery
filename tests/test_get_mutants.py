from qhery.get_mutants import mutantFinder
import pytest

working_dir='data'
sample_name='example'

@pytest.fixture()
def mf():
    return mutantFinder(working_dir=working_dir, sample_name='example')

def test_mutantFinder_init(mf):
    assert mf.bcftools_vcf == f'{working_dir}/{sample_name}.bcftools.vcf'
    assert mf.csq_file == f'{working_dir}/{sample_name}.csq.tsv'

@pytest.mark.parametrize(
    "working_dir,sample_name,expected",
    [
        ('data', 'example', ['E:T9I', 'M:A63T', 'M:Q19E', 'ORF3a:S165F', 'ORF3a:T223I', 'PLpro:G489S', 'PLpro:T24I', 'RdRP:P323L', 'S:D405N', 'S:D614G', 'S:D796Y', 'S:E484A', 'S:G142D', 'S:G339D', 'S:H655Y', 'S:I68T', 'S:K417N', 'S:L24S', 'S:N440K', 'S:N501Y', 'S:N679K', 'S:N764K', 'S:N969K', 'S:P681H', 'S:PPA25-27∆', 'S:Q493R', 'S:Q498R', 'S:Q954H', 'S:R408S', 'S:S371F', 'S:S373P', 'S:S375F', 'S:S477N', 'S:T19I', 'S:T376A', 'S:T478K', 'S:V213G', 'S:Y505H', '_3CLpro:P132H', 'nsp13:R392C', 'nsp1:S135R', 'nsp4:L264F', 'nsp4:L438F', 'nsp4:T327I', 'nsp4:T492I', 'nsp6:SGF106-108∆']),
        ('data', 'example2', ['E:T9I', 'M:A63T', 'M:Q19E', 'ORF3a:S165F', 'ORF3a:T223I', 'PLpro:G489S', 'PLpro:T24I', 'RdRP:P323L', 'S:D405N', 'S:D614G', 'S:D796Y', 'S:E484A', 'S:G142D', 'S:G339D', 'S:H655Y', 'S:I68T', 'S:K417N', 'S:L24S', 'S:N440K', 'S:N501Y', 'S:N679K', 'S:N764K', 'S:N969K', 'S:P681H', 'S:PPA25-27∆', 'S:Q493R', 'S:Q498R', 'S:Q954H', 'S:R408S', 'S:S371F', 'S:S373P', 'S:S375F', 'S:S477N', 'S:T19I', 'S:T376A', 'S:T478K', 'S:V213G', 'S:Y505H', '_3CLpro:P132H', 'nsp13:R392C', 'nsp1:S135R', 'nsp4:L264F', 'nsp4:L438F', 'nsp4:T327I', 'nsp4:T492I', 'nsp6:SGF106-108∆']),
        ('data', 'example3', ['E:T9I', 'M:A63T', 'M:Q19E', 'ORF3a:T223I', 'ORF6:D61L', 'ORF7b:L11F', 'PLpro:G489S', 'PLpro:S1675G', 'PLpro:T24I', 'RdRP:P323L', 'S:D405N', 'S:D614G', 'S:D796Y', 'S:E484A', 'S:F486V', 'S:G142D', 'S:G339D', 'S:H655Y', 'S:HV69-70∆', 'S:K417N', 'S:L24S', 'S:L452R', 'S:N440K', 'S:N501Y', 'S:N658S', 'S:N679K', 'S:N764K', 'S:N969K', 'S:P681H', 'S:PPA25-27∆', 'S:Q498R', 'S:Q954H', 'S:R408S', 'S:S371F', 'S:S373P', 'S:S375F', 'S:S477N', 'S:T19I', 'S:T376A', 'S:T478K', 'S:V213G', 'S:Y505H', '_3CLpro:P132H', 'nsp13:R392C', 'nsp14:I42V', 'nsp15:T112I', 'nsp1:KSF141-143∆', 'nsp1:S135R', 'nsp4:L264F', 'nsp4:T327I', 'nsp4:T492I', 'nsp6:SGF106-108∆']),
        ('data', 'example4', ['E:T9I', 'M:A63T', 'M:Q19E', 'ORF3a:T223I', 'ORF6:D61L', 'PLpro:G489S', 'PLpro:P192S', 'PLpro:T24I', 'RdRP:P323L', 'S:D405N', 'S:D614G', 'S:D796Y', 'S:E484A', 'S:G142D', 'S:G339D', 'S:H655Y', 'S:K417N', 'S:L24S', 'S:N440K', 'S:N501Y', 'S:N679K', 'S:N764K', 'S:N969K', 'S:P681H', 'S:PPA25-27∆', 'S:Q493R', 'S:Q498R', 'S:Q954H', 'S:R408S', 'S:S371F', 'S:S373P', 'S:S375F', 'S:S477N', 'S:T19I', 'S:T376A', 'S:T478K', 'S:V213G', 'S:Y505H', '_3CLpro:P132H', 'nsp13:R392C', 'nsp14:H26Q', 'nsp14:I42V', 'nsp1:S135R', 'nsp4:L264F', 'nsp4:L438F', 'nsp4:T327I', 'nsp4:T492I', 'nsp6:SGF106-108∆']),
        ('data', 'example5', ['E:T9I', 'M:A63T', 'M:Q19E', 'ORF3a:L140F', 'ORF3a:T223I', 'ORF3a:W131C', 'ORF6:D61L', 'ORF7a:E121D', 'PLpro:G489S', 'PLpro:T24I', 'RdRP:P323L', 'S:D405N', 'S:D614G', 'S:D796Y', 'S:E484A', 'S:G142D', 'S:G339D', 'S:H655Y', 'S:K417N', 'S:L24S', 'S:N440K', 'S:N501Y', 'S:N679K', 'S:N764K', 'S:N969K', 'S:P681H', 'S:PPA25-27∆', 'S:Q493R', 'S:Q498R', 'S:Q954H', 'S:R408S', 'S:S371F', 'S:S373P', 'S:S375F', 'S:S477N', 'S:T19I', 'S:T376A', 'S:T478K', 'S:V213G', 'S:Y505H', '_3CLpro:P132H', 'nsp13:R392C', 'nsp14:I42V', 'nsp15:T112I', 'nsp1:S135R', 'nsp4:A146V', 'nsp4:L264F', 'nsp4:L438F', 'nsp4:T327I', 'nsp4:T492I', 'nsp6:SGF106-108∆']),
    ],
)
def test_parse_csq(working_dir, sample_name,expected):
    mf = mutantFinder(working_dir=working_dir, sample_name=sample_name)
    res = mf.parse_csq()
    print(sample_name, sorted(res))
    assert all(a == b for a,b in zip(sorted(res), expected))