list_fn = ARGV[0]
list = File.open(list_fn,"r")
list.each do |line|
  endo = line.strip
  Dir.chdir(endo)
  if(File.exists?("seq_aa.hhm.lock"))
    system("rm -r seq_aa.hhm.lock")
  end
  if(File.exists?("seq_aa.hhm.log"))
    system("rm -r seq_aa.hhm.log")
  end
  system("~/src/cm_scripts/bin/make-hhm.pl  seq_aa.fasta -use_hhblits -n_procs 6 -diff 1000")
  system("~/src/hhsuite/hhsuite-2.0.2/bin/hhsearch  -cpu 6 -i seq_aa.hhm -d /work/brunette/src/endo/pdb_database/homing_endo_db")
  Dir.chdir("../..")
end


