################################################################
##         Bowtie2-Mapping-Single-ended.cwl                   ##
################################################################

cwlVersion: v1.0
class: CommandLineTool
doc: "[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Fast and sensitive read alignment"
label: Last updated 09/2019
hints:
  SoftwareRequirement:
    packages:
    - package: bowtie2
      version: [ 2.3.4.3 ]

################################################################
##           baseCommand and arguments definitions            ##
################################################################

baseCommand: [ bowtie2 ]

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.results_path) ]

arguments:
  - prefix: -S
    valueFrom: $(inputs.results_path.path)/$(inputs.SampleName).sam
  - prefix: -x
    valueFrom: $(inputs.idx_basedir.path)/$(inputs.idx_basename)
  - prefix: -k
    valueFrom: '50'
  - prefix: 
    valueFrom: --non-deterministic

################################################################
##               Inputs and Outputs Settings                  ##
################################################################

inputs:
  idx_basedir:
    label: "Path to the directory containing the index for the reference genome"
    type: Directory
  idx_basename:
    label: "Basename of the bowtie2 index files"
    type: string
  fq1:
    type: File
    label: "Comma-separated list of files containing unpaired reads to be aligned"
    inputBinding:
      prefix: -U
      itemSeparator: ","
  thread:
    label: "Launch NTHREADS parallel search threads"
    type: int
    inputBinding:
      prefix: -p
  SampleName:
    label: "Filename to write output to"
    type: string
  results_path:
    label: "Path to the results directory"
    type: Directory

outputs:
  bowtie2_sam:
    type: File
    outputBinding:
      glob: $(inputs.results_path.path)/$(inputs.SampleName).sam
