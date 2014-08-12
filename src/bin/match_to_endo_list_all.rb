lists = Dir.glob("list*")
lists.each do |list_fn|
    list_tmp = File.basename(list_fn)
    if(!File.exists?("w_#{list_tmp}"))
        system("touch w_#{list_tmp}")
        print("working on #{list_tmp}\n")
        list = File.open(list_fn,"r")
        list.each do |line|
            endo = line.strip
            if(File.exists?(endo))
                Dir.chdir(endo) 
                if(!File.exists?("seq_aa.hhr"))
                    if(File.exists?("seq_aa.hhm.lock"))
                        system("rm -r seq_aa.hhm.lock")
                    end
                    if(File.exists?("seq_aa.hhm.log"))
                        system("rm -r seq_aa.hhm.log")
                    end
                    system("~/src/cm_scripts/bin/make-hhm-fast.pl  seq_aa.fasta -use_hhblits -n_procs 6")
                    system("/work/robetta/src/rosetta_server/src/hhsuite-2.0.16-linux-x86_64/bin/hhsearch -cpu 6 -i seq_aa.hhm -d /work/brunette/src/endo/pdb_database/homing_endo_db")
                end
                Dir.chdir("../..")
            end
        end
    end
end
