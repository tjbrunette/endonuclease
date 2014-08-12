#endos = Dir.glob("endo_nonAmbig/*")
#endos = Dir.glob("endo_insufficient_DNA/*")
endos = Dir.glob("endo_ambig/*")
new_target_length = 40
current_target_length=15
endos.each do |endo|
  Dir.chdir(endo)
  upSeq_dna_f = File.open("upSeq_dna.fasta","r")
  upSeq_dna_target_f = File.open("closest_upSeq_target_dna.fasta","r")
  upSeq_dna_target_long_f = File.open("closest_upSeq_target_dna_40res.fasta","w")
  line=upSeq_dna_f.gets()
  upSeq_dna = upSeq_dna_f.gets().strip()
  upSeq_dna_target_long_f << upSeq_dna_target_f.gets()
  upSeq_dna_target =upSeq_dna_target_f.gets().strip()
  location = upSeq_dna.index(upSeq_dna_target)
  if(location == nil)
    print("upseq location doesn't exist for #{endo}\n")
  else
    upSeq_dna_target_longer= upSeq_dna[location-new_target_length+current_target_length..location+current_target_length]
    upSeq_dna_target_long_f << "#{upSeq_dna_target_longer}\n"
  end
  downSeq_dna_f = File.open("downSeq_dna.fasta","r")
  downSeq_dna_target_f = File.open("closest_downSeq_target_dna.fasta","r")
  downSeq_dna_target_long_f = File.open("closest_downSeq_target_dna_40res.fasta","w")
  line=downSeq_dna_f.gets()
  downSeq_dna = downSeq_dna_f.gets().strip()
  downSeq_dna_target_long_f << downSeq_dna_target_f.gets()
  downSeq_dna_target =downSeq_dna_target_f.gets().strip()
  location = downSeq_dna.index(downSeq_dna_target)
  if(location == nil)
    print("downseq location doesn't exist for #{endo}\n")
  else
    downSeq_dna_target_longer= downSeq_dna[location..location+new_target_length]
    downSeq_dna_target_long_f << "#{downSeq_dna_target_longer}\n"
  end
  upSeq_dna_f.close()
  upSeq_dna_target_f.close()
  upSeq_dna_target_long_f.close()
  downSeq_dna_f.close()
  downSeq_dna_target_f.close()
  downSeq_dna_target_long_f.close()
  Dir.chdir("../..")
end
