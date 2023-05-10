cwlVersion: v1.0
class: CommandLineTool
arguments:
- prefix: -hash_index
- prefix: -parse_seqids
baseCommand: makeblastdb
inputs:
  input_file:
    type: File
    inputBinding:
      prefix: -in
  database_name:
    type: File
    inputBinding:
      prefix: -out
  molecule_type:
    type: string
    inputBinding:
      prefix: -dbtype
outputs:
  out.pdb:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pdb
  out.phd:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).phd
  out.phi:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).phi
  out.phr:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).phr
  out.pin:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pin
  out.pog:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pog
  out.pos:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pos
  out.pot:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pot
  out.psq:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).psq
  out.ptf:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).ptf
  out.pto:
    type: File
    outputBinding:
      glob: $(inputs.results_path)/$(inputs.database_name.basename).pto
