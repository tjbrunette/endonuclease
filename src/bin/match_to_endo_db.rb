endos = Dir.glob("endo_nonAmbig/*")
endos += Dir.glob("endo_insufficient_DNA/*")
endos += Dir.glob("endo_ambig/*")
endos.each do |endo|
  Dir.chdir(endo)
  print("working on #{endo}\n")
  system("~/src/cm_scripts/bin/make-hhm.pl  seq_aa.fasta -use_hhblits -n_procs 20 -diff 1000")
  system("~/src/hhsuite/hhsuite-2.0.2/bin/hhsearch  -cpu 20 -i  seq_aa.hhm -d /work/brunette/src/endo/pdb_database/homing_endo_db")
  Dir.chdir("../..")
end


