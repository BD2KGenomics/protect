from __future__ import print_function
import logging
from pipelineWrapper import PipelineWrapperBuilder

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

desc = """UCSC Precision Immuno pipeline"""

config = ("""patients:
{samples}

Universal_Options:
    dockerhub: aarjunrao
    java_Xmx: 20G
    sse_key: {sse_key}
    sse_key_is_master: {sse_key_is_master}
    cghub_key:
    storage_location: {storage-location}
    output_folder: {output-folder}

# These options are for each specific tool. You probably don't need to change any of this!
cutadapt:
    a : AGATCGGAAGAG
    A : AGATCGGAAGAG

star:
    type : star # use starlong if your reads are > 150bp
    tool_index : {star_index}

bwa:
    tool_index : {bwa_index}

rsem:
    tool_index : {rsem_index}

mut_callers:
    genome_fasta : {genome_fasta}
    genome_fai : {genome_fai}
    genome_dict : {genome_dict}
    cosmic_vcf : {cosmic_vcf}
    cosmic_idx : {cosmic_idx}
    dbsnp_vcf : {dbsnp_vcf}
    dbsnp_idx : {dbsnp_idx}
    dbsnp_tbi : {dbsnp_tbi}
    java_Xmx : 2G
    strelka_config: {strelka_config}

snpeff:
    tool_index : {snpeff}
    java_Xmx : 20G

transgene:
    gencode_peptide_fasta : {transgene}

phlat:
    tool_index : {phlat}

mhci:
    method_file : {mhci}
    pred : IEDB_recommended

mhcii:
    method_file : {mhcii}
    pred : IEDB_recommended

rank_boost:
    mhci_combo : 0.0,0.0,0.32,0.0,0.68
    mhcii_combo : 0.25,0.25,0.25,0.25

mhc_pathway_assessment:
    genes_file : {mhc_pathway_assessment}""")

if __name__ == '__main__':
    wrapper = PipelineWrapperBuilder('ProTECT', desc, config)
    parser = wrapper.get_args()

    parser.add_argument('--tumor-dna', default=[], action="append",
                        help='Path for the tumor fastq.')
    parser.add_argument('--normal-dna',  default=[], action="append",
                        help='Path for the normal fastq.')
    parser.add_argument('--tumor-rna',  default=[], action="append",
                        help='Path for the tumor RNA fastq.')
    parser.add_argument('--tumor-dna2',  default=[], action="append",
                        help='Path for the tumor fastq pair.')
    parser.add_argument('--normal-dna2',  default=[], action="append",
                        help='Path for the normal fastq.')
    parser.add_argument('--tumor-rna2',  default=[], action="append",
                        help='Path for the tumor RNA fastq.')

    parser.add_argument('--star-index', type=str, default="S3://cgl-protect-data/hg19_references/star_100_indexes.tar.gz",
                        help='Use starlong if your reads are > 150bp')

    parser.add_argument('--bwa-index', type=str, default="S3://cgl-protect-data/hg19_references/bwa_index.tar.gz")

    parser.add_argument('--rsem-index', type=str, default="S3://cgl-protect-data/hg19_references/rsem_index.tar.gz")

    parser.add_argument('--genome-fasta', type=str, default="S3://cgl-protect-data/hg19_references/hg19.fa.tar.gz")
    parser.add_argument('--genome-fai', type=str, default="S3://cgl-protect-data/hg19_references/hg19.fa.fai.tar.gz")
    parser.add_argument('--genome-dict', type=str, default="S3://cgl-protect-data/hg19_references/hg19.dict.tar.gz")

    parser.add_argument('--cosmic-vcf', type=str, default="S3://cgl-protect-data/hg19_references/CosmicCodingMuts.vcf.tar.gz")
    parser.add_argument('--cosmic-idx', type=str, default="S3://cgl-protect-data/hg19_references/CosmicCodingMuts.vcf.idx.tar.gz")

    parser.add_argument('--dbsnp-vcf', type=str, default="S3://cgl-protect-data/hg19_references/dbsnp_coding.vcf.gz")
    parser.add_argument('--dbsnp-idx', type=str, default="S3://cgl-protect-data/hg19_references/dbsnp_coding.vcf.idx.tar.gz")
    parser.add_argument('--dbsnp-tbi', type=str, default="S3://cgl-protect-data/hg19_references/dbsnp_coding.vcf.gz.tbi")

    parser.add_argument('--strelka-config', type=str, default="S3://cgl-protect-data/hg19_references/strelka_bwa_WXS_config.ini.tar.gz")

    parser.add_argument('--snpeff', type=str, default="S3://cgl-protect-data/hg19_references/snpeff_index.tar.gz")

    parser.add_argument('--transgene', type=str, default="S3://cgl-protect-data/hg19_references/gencode.v19.pc_translations_NOPARY.fa.tar.gz")

    parser.add_argument('--phlat', type=str, default="S3://cgl-protect-data/hg19_references/phlat_index.tar.gz")

    parser.add_argument('--mhci', type=str, default="S3://cgl-protect-data/hg19_references/mhci_restrictions.json.tar.gz")
    parser.add_argument('--mhcii', type=str, default="S3://cgl-protect-data/hg19_references/mhcii_restrictions.json.tar.gz")

    parser.add_argument('--mhc-pathway-assessment', type=str, default="S3://cgl-protect-data/hg19_references/mhc_pathway_genes.json.tar.gz")

    parser.add_argument('--sse-key-is-master', type=str, default='False',
                        help='Indicates if the passed sse-key is the master key.')

    parser.add_argument('--sse-key', type=str, default='',
                        help='Path to the desired SSE-key, if any.')

    parser.add_argument('--autoscale', dest='autoscale', default=False,
                        help="Indicates whether to use Toil's autoscaling capabilities")

    parser.add_argument('--work-mount', required=True,
                        help='Mount where intermediate files should be written. This directory '
                             'should be mirror mounted into the container.')

    parser.add_argument('--provisioner', type=str, default='aws',
                        help='Sets the location to where the autoscale should be held')

    parser.add_argument('--nodeType', type=str, default='t2.micro',
                        help="sets the type of node used to help with the clustering")

    args = parser.parse_args()
    command = []
    wrapper.run(args, command)
