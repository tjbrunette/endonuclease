#endos = Dir.glob("endo_nonAmbig/*")
#endos += Dir.glob("endo_insufficient_DNA/*")
#endos += Dir.glob("endo_ambig/*")
endos = ["endo_nonAmbig/AJ841808.1_2530.CAH56513.1","endo_nonAmbig/AJ841805.1_2483.CAH56510.1"]
endos.each do |endo|
  Dir.chdir(endo)
  print("working on #{endo}\n")
  system("~/src/endo/cm_scripts/bin/convert_aln.pl seq_aa.hhr -format_in hhsearch -format_out grishin -max_templates 1 > seq_aa.filt")
  system("~/src/rosetta/rosetta_source/bin/partial_thread.default.linuxgccrelease -in:file:alignment seq_aa.filt -cm:aln_format grishin -in:file:fasta seq_aa.fasta -in:file:template_pdb /work/brunette/src/endo/pdb_database/template_proteins_no_dna/*.pdb")
  Dir.chdir("../..")
end


