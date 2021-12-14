cwlVersion: v1.2
class: Workflow
label: Iterative LD pruning by group
doc: This workflow does LD pruning iteratively by group, keeping only variants not in LD in any group.
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
  - id: strata_file
    label: strata file
    doc: RDS file with list of sample.id
    type: File
    sbg:fileTypes: RDS
  - id: variant_file
    label: variant file
    doc: RDS file with list of variant.id
    type: File
    sbg:fileTypes: RDS
  - id: maf
    label: MAF threshold
    doc: MAF threshold for pruning in pooled sample
    type: float?
    sbg:toolDefaultValue: '0.01'
  - id: ld_threshold
    label: LD threshold
    doc: LD pruning |r| threshold
    type: float?
    sbg:toolDefaultValue: sqrt(0.1)
  - id: out_prefix
    label: output prefix
    doc: output prefix (will have ".rds" appended)
    type: string?
    default: var_eff
    sbg:toolDefaultValue: pruned

outputs:
  - id: output_file
    label: Output file
    doc: RDS file containing vector of variant.id
    type: File
    outputSource: [ cwl_tool_wrapper_workflow_step/output_file ]

steps:
  cwl_tool_wrapper_workflow_step:
    in:
      gds_file: gds_file
      strata_file: strata_file
      variant_file: variant_file
      maf: maf
      ld_threshold: ld_threshold
      out_prefix: out_prefix
      
    out:
      - id: output_file

    run: 
      cwlVersion: v1.2
      class: CommandLineTool
      label: Iterative LD pruning by group
      doc: This tool does LD pruning iteratively by group, keeping only variants not in LD in any group.

      requirements:
      - class: ShellCommandRequirement
      - class: DockerRequirement
        dockerPull: uwgac/simphen:0.2.2
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
        sbg:fileTypes: GDS
      - id: strata_file
        label: strata file
        doc: RDS file with list of sample.id
        type: File
        inputBinding:
          prefix: --strata_file
          position: 2
          shellQuote: false
        sbg:fileTypes: RDS
      - id: variant_file
        label: variant file
        doc: RDS file with list of variant.id
        type: File
        inputBinding:
          prefix: --variant_file
          position: 3
          shellQuote: false
        sbg:fileTypes: RDS
      - id: maf
        label: MAF threshold
        doc: MAF threshold for pruning in pooled sample
        type: float?
        inputBinding:
          prefix: --maf
          position: 4
          shellQuote: false
        sbg:toolDefaultValue: '0.01'
      - id: ld_threshold
        label: LD threshold
        doc: LD pruning |r| threshold
        type: float?
        inputBinding:
          prefix: --ld_threshold
          position: 5
          shellQuote: false
        sbg:toolDefaultValue: sqrt(0.1)
      - id: out_prefix
        label: output prefix
        doc: output prefix (will have ".rds" appended)
        type: string?
        default: var_eff
        inputBinding:
          prefix: --out_file
          position: 6
          valueFrom: |-
            ${
                var chr = inputs.gds_file.nameroot.split('chr')[1]
                return inputs.out_prefix + "_pruned_chr" + chr + ".rds"
            }
          shellQuote: false
        sbg:toolDefaultValue: pruned

      outputs:
      - id: output_file
        label: Output file
        doc: RDS file containing vector of variant.id
        type: File
        outputBinding:
          glob: '*.rds'

      baseCommand:
      - Rscript
      - /usr/local/simulate_phenotypes/tools/iterative_ld_pruning.R

      stdout: job.out.log
      hints:
      - class: sbg:SaveLogs
        value: job.out.log
