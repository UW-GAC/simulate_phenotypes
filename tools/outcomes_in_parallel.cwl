cwlVersion: v1.0
class: CommandLineTool
label: Generate random phenotypes from covariance blocks
doc: This tool generates random correlated phenotypes based on a block diagonal covariance matrix. Block indices should be computed in advance using the function simphen:::block_indices, and the covariance matrix and block indices should be saved in a single RData object. See https://github.com/UW-GAC/simulate_phenotypes/blob/master/notebooks/generate_outcomes.Rmd for an example.
$namespaces:
  sbg: https://sevenbridges.com

requirements:
- class: ShellCommandRequirement
- class: ResourceRequirement
  coresMin: 0
  ramMin: 4000
- class: DockerRequirement
  dockerPull: uwgac/simphen:0.2.2
- class: InlineJavascriptRequirement

inputs:
- id: covmat_file
  label: covmat file
  doc: RData file with objects covmat and blocks
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
    valueFrom: ${ return inputs.out_prefix + ".rds"}
    shellQuote: false
- id: varComp1
  label: variance component 1
  doc: genetic variant component
  type: int
  inputBinding:
    prefix: --varComp1
    position: 3
    shellQuote: false
- id: varComp2
  label: variance component 2
  doc: error variant component
  type: int
  inputBinding:
    prefix: --varComp2
    position: 4
    shellQuote: false
- id: num_outcomes
  label: number of outcomes
  doc: number of outcomes to generate
  type: int
  inputBinding:
    prefix: --num_outcomes
    position: 5
    shellQuote: false

outputs:
- id: output_file
  type: File
  outputBinding:
    glob: '*.rds'

baseCommand:
- Rscript
- /usr/local/simulate_phenotypes/tools/outcomes_in_parallel.R

stdout: job.out.log
hints:
- class: sbg:SaveLogs
  value: job.out.log
