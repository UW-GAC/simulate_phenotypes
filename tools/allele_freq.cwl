cwlVersion: v1.1
class: CommandLineTool
label: Allele frequency
doc: This tool calculates allele frequency for variants in a GDS file.
$namespaces:
  sbg: https://sevenbridges.com
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/

dct:creator:
  "@id": "https://orcid.org/0000-0002-7231-9745"
  foaf:name: Stephanie Gogarten
  foaf:mbox: "mailto:sdmorris@uw.edu"

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: ${ return inputs.cpu }
- class: DockerRequirement
  dockerPull: uwgac/simphen:0.2.2
- class: InlineJavascriptRequirement

inputs:
- id: gds_file
  label: GDS file
  doc: GDS file
  type: File
  inputBinding:
    position: 1
    shellQuote: false
- id: out_prefix
  label: Output prefix
  doc: Prefix for output files
  type: string
  inputBinding:
    position: 2
    valueFrom: |-
      ${
          var chr = inputs.gds_file.nameroot.split('chr')[1]
          return inputs.out_prefix + "_allele_freq_chr" + chr + ".rds"
      }
    shellQuote: false
- id: sample_file
  label: Sample file
  doc: RDS file with vector of sample.id to include in calculating allele frequency
  type: File?
  inputBinding:
    prefix: --sample.file
    position: 3
    shellQuote: false
- id: cpu
  type: int?
  inputBinding:
    prefix: --cpu
    position: 4
    shellQuote: false
  sbg:toolDefaultValue: '1'

outputs:
- id: output_file
  label: Output file
  doc: |-
    RDS file containing R data.frame with variant.id, chromosome, position, and alternate allele frequency.
  type: File
  outputBinding:
    glob: '*.rds'
stdout: job.out.log

baseCommand:
- Rscript
- /usr/local/simulate_phenotypes/tools/allele_freq.R

hints:
- class: sbg:SaveLogs
  value: job.out.log
