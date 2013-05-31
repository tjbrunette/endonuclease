endos = Dir.glob("endo_nonAmbig/*")
endos += Dir.glob("endo_insufficient_DNA/*")
endos += Dir.glob("endo_ambig/*")
Dir.mkdir("duplicates")
Dir.mkdir("duplicates/endo_nonAmbig")
Dir.mkdir("duplicates/endo_insufficient_DNA")
Dir.mkdir("duplicates/endo_ambig")
aaSeq_name_h = Hash.new
endos.each do |endo|
  fasta = File.open("#{endo}/seq_aa.fasta")
  name = fasta.gets.lstrip
  aaSeq = fasta.gets.lstrip
  fasta.close
  if(aaSeq_name_h.has_key?(aaSeq))
    real_name = endo.split("/")[1]
    type = endo.split("/")[0]
    print("duplicates:#{name} and #{aaSeq_name_h[aaSeq]}\n")
    system("mv #{endo} duplicates/#{type}/#{real_name}")
  else
    aaSeq_name_h[aaSeq]=name
  end
end


