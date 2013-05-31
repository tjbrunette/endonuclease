def convert_string_to_regex(aa_seq)
  aa_regex = ""
  aa_seq.each_char do |t|
    aa_regex += "#{t}(-)*"
  end
  return(Regexp.new(aa_regex))
end

def get_max_length_each_sheet(template_endos)
  sheet_length_a = Array.new
  template_endos.each do |key,aa_strings|
    for ii in 0..aa_strings.size-1
      if(sheet_length_a.size <= ii)
        sheet_length_a.push(aa_strings[ii].size)
      else
        if(sheet_length_a[ii]<aa_strings[ii].size)
          sheet_length_a[ii] = aa_strings[ii].size
        end
      end
    end
  end
  return(sheet_length_a)
end

def remove_dashes_from_match(target,start_res,stop_res,match)
  final_match = ""
  for ii in 0..stop_res-start_res
    if(match[ii] != "-")
      final_match += target[start_res+ii]
    end
  end
  return(final_match)
end

def pad_dashes(current_match,desired_length)
  final_match = current_match
  for ii in current_match.size..desired_length-1
    final_match += "-"
  end
  return(final_match)
end

def get_beta_residues_equal_spacing(target,template,template_endos,pdb,sheet_length_a)
  beta_residues = ""
  print("#{pdb}\n")
  sheet_index = 0
  template_endos[pdb].each do |seq|
    regexp_aa = convert_string_to_regex(seq)
    template_parts = template.partition(regexp_aa)
    start_res = template.index(regexp_aa)
    if(start_res != nil)
      match = template_parts[1]
      stop_res = start_res+match.size-1
      match_no_dashes = remove_dashes_from_match(target,start_res,stop_res,match)
      match_correct_length = pad_dashes(match_no_dashes,sheet_length_a[sheet_index])
      beta_residues += "#{match_correct_length}"
    else
      match_correct_length = pad_dashes("",sheet_length_a[sheet_index])
      beta_residues += "#{match_correct_length}"
    end
    sheet_index += 1
  end
  return(beta_residues);
end


endos = Dir.glob("endo_nonAmbig/*")
endos += Dir.glob("endo_insufficient_DNA/*")
endos += Dir.glob("endo_ambig/*")
beta_sheets_out = File.open("beta_sheets.out","w")
#endos = ["endo_nonAmbig/AJ841808.1_2530.CAH56513.1","endo_nonAmbig/AJ841805.1_2483.CAH56510.1"]
#template_endos = Hash["3ool"=>["YIR","CMQFEW","HKKERV","LVITWGAQ","GGK","IVL","CYVKIN","KPIIYI"], "3qqy"=>["SFLLRIRN","YSTELGFQITLH","VIAN","AVSLKV","GCFFVNLIK","VQVQLVFSITQ","YIKEKNK","FSWLDFVV"]]
template_endos = Hash["2qoj"=>["GYFSITKK","YLTYELGIELS","IVSFR","MVALRI","GCFS","SFDIAQR","YL","CSKLKV"],"2vbj"=>["GSIIAQIEP","HRLKLTFKVTQK","YVRD","VSNYIL","","","",""],"3e54"=>["GSFSVSIKF","VRLDPVFSITQ","RIME","YVYVV","","","",""],"3fd2"=>["GSIYAKLIP","YQVSLAISFIQR","NLR","IADYTII","SIYAKLIPR","QVSLAISFIQRK","LRK","ADYTIIG"],"3ool"=>["YIR","CMQFEW","HKKERV","LVITWGAQ","GGK","IVL","CYVKIN","KPIIYI"],"3qqy"=>["SFLLRIRN","YSTELGFQITLH","VIAN","AVSLKV","CFFVNLIKS","QVQLVFSITQH","IKEKNKS","SWLDFVVT"],"3r7p"=>["SFMLTVSK","WSVRPRFRIGL","IITS","ARIRF","GSFYIRIAK","YQVQSVFQITQD","NIRIR","VDLVVT"]]
sheet_length_a = get_max_length_each_sheet(template_endos)
endos.each do |endo|
  Dir.chdir(endo)
  print("working on #{endo}\n")
  system("~/src/endo/cm_scripts/bin/convert_aln.pl seq_aa.hhr -format_in hhsearch -format_out grishin -max_templates 100 -nosort > seq_aa.filt")
  endo_filt = File.open("seq_aa.filt","r")
  done = false
  beta_residues = ""
  while(!done && line = endo_filt.gets)
    position = line.split(" ")
    if(position[0] == "##")
      pdb = position[2][0..3]
      if(template_endos.include?(pdb))
        done = true
        line=endo_filt.gets
        line=endo_filt.gets
        line=endo_filt.gets
        target = line.split(" ")[1]
        line = endo_filt.gets
        template = line.split(" ")[1]
        beta_residues = get_beta_residues_equal_spacing(target,template,template_endos,pdb,sheet_length_a)
      end
    end
  end
  beta_sheets_out << ">#{endo}\n"
  beta_sheets_out << "#{beta_residues}\n"
  Dir.chdir("../..")
end
