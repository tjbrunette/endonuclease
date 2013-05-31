#convert fasta to pssm file
def construct_pssm(fasta_file_name)
  fasta_file = File.new(fasta_file_name,'r')
  pssm_file = File.new("#{fasta_file_name}.pssm",'w')
  line = fasta_file.gets()
  line = fasta_file.gets()
  pssm_file << "key  a    c    g    t\n"
  for i in 0..(line.size()-2)
    pssm_hash = Hash.new
    pssm_hash['a'] = 0.0
    pssm_hash['c'] = 0.0
    pssm_hash['g'] = 0.0
    pssm_hash['t'] = 0.0
    if ((i<3)||(i>11))
      pssm_hash[line[i..i]] = -1.0
    else
      pssm_hash[line[i..i]] = -1.0
    end
    pssm_file << "#{i}  #{pssm_hash['a']}  #{pssm_hash['c']}  #{pssm_hash['g']}  #{pssm_hash['t']}\n"
  end
  pssm_file.close()
end

def identify_target_site(neighbor_target_site_fn,flanking_res_fn,pred_target_site_fn)
  construct_pssm(neighbor_target_site_fn)
  system("~/src/endo/bin/pssm -s #{flanking_res_fn} -p #{neighbor_target_site_fn}.pssm > #{pred_target_site_fn}")
end

#extracts files from predicted pssm file output and puts it in the output
def extract_pred_target_DNA(pred_file_name,targetDna_file_name)
  pred_file = File.new(pred_file_name,'r')
  last_seq = ""
  proteinName = ""
  while(line = pred_file.gets())
    position = line.split(" ")
    if(position.size()>2)
      if(position[2][0..3] ==">ref")
        proteinName = position[2]
        last_seq = position[1]
      end
    end
  end
  targetDna_file = File.new(targetDna_file_name,"w")
  targetDna_file << "#{proteinName}\n"
  targetDna_file << "#{last_seq}\n"
  targetDna_file.close()
end

Dir.chdir("endo_ambig")
endos = Dir.glob("*")
endos.each do |f|
  Dir.chdir(f)
  identify_target_site("closest_downSeq_target_dna.fasta","downSeq_dna.fasta","downSeq_target_dna_pred")
  identify_target_site("closest_upSeq_target_dna.fasta","upSeq_dna.fasta","upSeq_target_dna_pred")
  Dir.chdir("..")
end
