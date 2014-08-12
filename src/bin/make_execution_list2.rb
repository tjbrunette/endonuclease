#endos = Dir.glob("endo_nonAmbig/*")
endos = Dir.glob("endo_insufficient_DNA/*")
#endos += Dir.glob("endo_ambig/*")
list2 = File.open("day4List.txt","w")
endos.each do |endo|
  if(!File.exists?("#{endo}/seq_aa.hhr"))
    list2 << "#{endo}\n"
  end
end
list2.close
