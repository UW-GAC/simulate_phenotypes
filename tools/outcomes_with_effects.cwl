cwlVersion: v1.2
class: Workflow
label: Simulate phenotypes with effects
doc: This workflow generates random correlated phenotypes based on a block diagonal covariance matrix, then simulates the genetic effects of correlated variants and adds them to the outcome.
$namespaces:
  sbg: https://sevenbridges.com
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  
dct:creator:
  "@id": "http://orcid.org/0000-0002-7231-9745"
  foaf:name: Stephanie Gogarten
  foaf:mbox: sdmorris@uw.edu

inputs:
- id: gds_file
  label: GDS file
  doc: GDS file
  type: File
  sbg:fileTypes: GDS
- id: covmat_file
  label: covmat file
  doc: RDS file with block-diagonal covariance matrix
  type: File
  sbg:fileTypes: RDS
- id: variant_file
  label: variant file
  doc: RDS file with variant.id, as either a vector or a data.frame with two columns, one called 'variant.id' and the other called 'h2' or 'beta'
  type: File
  sbg:fileTypes: RDS
- id: h2
  label: heritability
  doc: heritability (used if variant.id is a vector)
  type: float?
- id: beta
  label: beta
  doc: beta (used if variant.id is a vector and h2 is not supplied)
  type: float?
- id: varComp1
  label: genetic variance component
  doc: genetic variance component
  type: int
- id: varComp2
  label: error variance component
  doc: error variance component
  type: int
- id: seed
  label: seed
  doc: seed for random outcome generation
  type: int?
- id: out_prefix
  label: output prefix
  doc: output prefix (will have ".rds" appended)
  type: string?
  sbg:toolDefaultValue: outcome_var_eff

outputs:
- id: output_file
  label: Output file
  doc: RDS file containing vector of outcomes
  type: File
  outputSource: [ cwl_tool_wrapper_workflow_step/output_file ]

steps:
  cwl_tool_wrapper_workflow_step:
    in:
      gds_file: gds_file
      covmat_file: covmat_file
      variant_file: variant_file
      h2: h2
      beta: beta
      varComp1: varComp1
      varComp2: varComp2
      seed: seed
      out_prefix: out_prefix

    out:
      - id: output_file

    run:
      cwlVersion: v1.2
      class: CommandLineTool
      label: Simulate phenotypes with effects
      doc: This workflow generates random correlated phenotypes based on a block diagonal covariance matrix, then simulates the genetic effects of correlated variants and adds them to the outcome.

      requirements:
      - class: ShellCommandRequirement
      - class: ResourceRequirement
        coresMin: 0
        ramMin: 4000
      - class: DockerRequirement
        dockerPull: uwgac/simphen:0.2.2-1
      - class: InlineJavascriptRequirement

      inputs:
      - id: gds_file
        label: GDS file
        doc: GDS file
        type: File
        inputBinding:
          prefix: --gds_file
          position: 1
          shellQuote: false
      - id: covmat_file
        label: covmat file
        doc: RDS file with block-diagonal covariance matrix
        type: File
        inputBinding:
          prefix: --covmat_file
          position: 2
          shellQuote: false
      - id: variant_file
        label: variant file
        doc: RDS file with variant.id, as either a vector or a data.frame with two columns, one called 'variant.id' and the other called 'h2' or 'beta'
        type: File
        inputBinding:
          prefix: --variant_file
          position: 3
          shellQuote: false
        sbg:fileTypes: RDS
      - id: h2
        label: heritability
        doc: heritability for each stratum
        type: float?
        inputBinding:
          prefix: --h2
          position: 4
          shellQuote: false
      - id: beta
        label: beta
        doc: beta for each stratum
        type: float?
        inputBinding:
          prefix: --beta
          position: 5
          shellQuote: false
      - id: varComp1
        label: genetic variance component
        doc: genetic variance component
        type: int
        inputBinding:
          prefix: --varComp1
          position: 6
          shellQuote: false
      - id: varComp2
        label: error variance component
        doc: error variance component
        type: int
        inputBinding:
          prefix: --varComp2
          position: 7
          shellQuote: false
      - id: seed
        label: seed
        doc: seed for random outcome generation
        type: int?
        inputBinding:
          prefix: --seed
          position: 8
          shellQuote: false
      - id: out_prefix
        label: output prefix
        doc: output prefix (will have ".rds" appended)
        type: string?
        default: outcome_var_eff
        sbg:toolDefaultValue: outcome_var_eff
        inputBinding:
          prefix: --out_file
          position: 9
          valueFrom: ${ return inputs.out_prefix + ".rds"}
          shellQuote: false

      outputs:
      - id: output_file
        label: Output file
        doc: RDS file containing vector of outcomes
        type: File
        outputBinding:
          glob: '*.rds'

      baseCommand:
      - Rscript
      - /usr/local/simulate_phenotypes/tools/outcomes_with_effects.R

      stdout: job.out.log
      hints:
      - class: sbg:SaveLogs
        value: job.out.log

