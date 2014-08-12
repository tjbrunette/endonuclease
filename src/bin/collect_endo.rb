endos = Dir.glob("endo_nonAmbig/*")
endos += Dir.glob("endo_ambig/*")
endos += Dir.glob("endo_insufficient_DNA/*")
#Dir.mkdir("non_endo")
#Dir.mkdir("non_endo/endo_nonAmbig")
#Dir.mkdir("non_endo/endo_insufficient_DNA")
#Dir.mkdir("non_endo/endo_ambig")
endos.each do |endo|
  type = endo.split("/")[0]
  real_name = endo.split("/")[1]
  if(File.exists?("#{endo}/seq_aa.hhr"))
    hhr = File.open("#{endo}/seq_aa.hhr")
    found_desired = false
    prob = -1
    while((line = hhr.gets) && (found_desired == false))
      position = line.split(" ")
      if(position[0] == "1")
        found_desired = true
        prob = position[2].to_f
      end
    end
    hhr.close
    if(prob >= 50)
      print("typeA #{prob} mv #{endo} endos/#{type}/#{real_name}\n")
      system("mv #{endo} endos/#{type}/#{real_name}")
    end
  end
end

