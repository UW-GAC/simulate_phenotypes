cwlVersion: v1.2
class: Workflow
label: Simulate variant effects
doc: This workflow simulates the genetic effects of correlated variants in groups of samples separately and in all samples pooled. The "simulate_correlated_outcomes" tool should be run first to simulate outcomes in the absence of genetic effects, and both the output of that tool and the input file with covariance matrix and block indices should be provided here as well. This tool will randomly select one of the previously simulated outcomes, and num_variants from the supplied list of variants. Variants may be tested in blocks, where the effects of multiple variants are added together. This saves compuatation time in running the association test, but it requires that the variants selected are not in LD in any of the groups being tested (use the "iterative_ld_pruning" tool first to select variants). See the vignette of the simphen package (https://github.com/UW-GAC/simulate_phenotypes/blob/master/simphen/vignettes/simulate_phenotypes.Rmd) for a description of this suite of tools. 
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
  doc: RData file with objects covmat and blocks
  type: File
  sbg:fileTypes: RDATA
- id: covar_file
  label: covar file
  doc: RDS file with AnnotatedDataFrame containing sample.id and covariates
  type: File
  sbg:fileTypes: RDS
- id: outcome_file
  label: outcome file
  doc: RDS file with random outcomes
  type: File
  sbg:fileTypes: RDS
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
- id: h2
  label: heritability
  doc: heritability for each stratum
  type: float[]?
- id: beta
  label: beta
  doc: beta for each stratum
  type: float[]?
- id: varComp1
  label: genetic variance component
  doc: genetic variance component
  type: int
- id: varComp2
  label: error variance component
  doc: error variance component
  type: int
- id: num_variants
  label: number of variants
  doc: number of variants to test
  type: int?
- id: variant_block_size
  label: variant block size
  doc: number of variants to add to each outcome at the same time
  type: int?
  sbg:toolDefaultValue: '10'
- id: out_prefix
  label: output prefix
  doc: output prefix (will have ".rds" appended)
  type: string?
  sbg:toolDefaultValue: var_eff
- id: num_cores
  type: int?
  sbg:toolDefaultValue: 1

outputs:
- id: output
  label: Output file
  doc: RDS file containing data.frame of test results
  type: File
  outputSource: [ cwl_tool_wrapper_workflow_step/output_file ]

steps:
  cwl_tool_wrapper_workflow_step:
    in:
      gds_file: gds_file
      covmat_file: covmat_file
      covar_file: covar_file
      outcome_file: outcome_file
      strata_file: strata_file
      variant_file: variant_file
      h2: h2
      beta: beta
      varComp1: varComp1
      varComp2: varComp2
      num_variants: num_variants
      variant_block_size: variant_block_size
      out_prefix: out_prefix
      num_cores: num_cores

    out:
      - id: output_file

    run:
      cwlVersion: v1.2
      class: CommandLineTool
      label: Simulate variant effects
      doc: This tool simulates the genetic effects of correlated variants in groups of samples separately and in all samples pooled. The "simulate_correlated_outcomes" tool should be run first to simulate outcomes in the absence of genetic effects, and both the output of that tool and the input file with covariance matrix and block indices should be provided here as well. This tool will randomly select one of the previously simulated outcomes, and num_variants from the supplied list of variants. Variants may be tested in blocks, where the effects of multiple variants are added together. This saves compuatation time in running the association test, but it requires that the variants selected are not in LD in any of the groups being tested (use the "iterative_ld_pruning" tool first to select variants). See the vignette of the simphen package (https://github.com/UW-GAC/simulate_phenotypes/blob/master/simphen/vignettes/simulate_phenotypes.Rmd) for a description of this suite of tools. 

      requirements:
      - class: ShellCommandRequirement
      - class: ResourceRequirement
        coresMin: ${ return inputs.num_cores }
        ramMin: 8000
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
      - id: covmat_file
        label: covmat file
        doc: RData file with objects covmat and blocks
        type: File
        inputBinding:
          prefix: --covmat_file
          position: 2
          shellQuote: false
        sbg:fileTypes: RDATA
      - id: covar_file
        label: covar file
        doc: RDS file with AnnotatedDataFrame containing sample.id and covariates
        type: File
        inputBinding:
          prefix: --covar_file
          position: 3
          shellQuote: false
        sbg:fileTypes: RDS
      - id: outcome_file
        label: outcome file
        doc: RDS file with random outcomes
        type: File
        inputBinding:
          prefix: --outcome_file
          position: 4
          shellQuote: false
        sbg:fileTypes: RDS
      - id: strata_file
        label: strata file
        doc: RDS file with list of sample.id
        type: File
        inputBinding:
          prefix: --strata_file
          position: 5
          shellQuote: false
        sbg:fileTypes: RDS
      - id: variant_file
        label: variant file
        doc: RDS file with list of variant.id
        type: File
        inputBinding:
          prefix: --variant_file
          position: 6
          shellQuote: false
        sbg:fileTypes: RDS
      - id: h2
        label: heritability
        doc: heritability for each stratum
        type: float[]?
        inputBinding:
          prefix: --h2
          position: 7
          shellQuote: false
      - id: beta
        label: beta
        doc: beta for each stratum
        type: float[]?
        inputBinding:
          prefix: --beta
          position: 8
          shellQuote: false
      - id: varComp1
        label: genetic variance component
        doc: genetic variance component
        type: int
        inputBinding:
          prefix: --varComp1
          position: 9
          shellQuote: false
      - id: varComp2
        label: error variance component
        doc: error variance component
        type: int
        inputBinding:
          prefix: --varComp2
          position: 10
          shellQuote: false
      - id: num_variants
        label: number of variants
        doc: number of variants to test
        type: int?
        inputBinding:
          prefix: --num_variants
          position: 11
          shellQuote: false
      - id: variant_block_size
        label: variant block size
        doc: number of variants to add to each outcome at the same time
        type: int?
        default: 10
        inputBinding:
          prefix: --variant_block_size
          position: 12
          shellQuote: false
        sbg:toolDefaultValue: '10'
      - id: out_prefix
        label: output prefix
        doc: output prefix (will have ".rds" appended)
        type: string?
        default: var_eff
        inputBinding:
          prefix: --out_file
          position: 13
          valueFrom: |-
            ${
                var chr = inputs.gds_file.nameroot.split('chr')[1]
                return inputs.out_prefix + "_effects_chr" + chr + ".rds"
            }
          shellQuote: false
        sbg:toolDefaultValue: var_eff
      - id: num_cores
        type: int?
        default: 1
        inputBinding:
          prefix: --num_cores
          position: 14
          shellQuote: false

      outputs:
      - id: output_file
        type: File
        outputBinding:
          glob: '*.rds'

      baseCommand:
      - Rscript
      - /usr/local/simulate_phenotypes/tools/test_variant_effects.R

      stdout: job.out.log
      hints:
      - class: sbg:SaveLogs
        value: job.out.log
