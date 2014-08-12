#!/usr/bin/ruby

#match_ambig_2_closestNonAmbig()
def match_ambig_2_closestNonAmbig(blast_dir,basepath)
  #construct nonAmbig blast db and store with nonAmbig
  #construct_nonAmbig_blastdb(blast_dir)
  blast_each_ambig(basepath)
  construct_pssm_each_nonAmbig()
  identify_target_each_ambig(basepath)
end

#run secondary_structure_prediction on all fasta files
def predict_secondary_structure_all()
  non_ambig_dir = "endo_nonAmbig"
  proteins=Dir.glob("#{non_ambig_dir}/*/seq_aa.fasta")
  fastas_to_combine = ""
  proteins.each do |f|
    fastas_to_combine = "#{fastas_to_combine} #{f}"
  end
  system ("cat #{fastas_to_combine} > #{non_ambig_dir}/endo_all_nonAmbig.fasta\n")
  system ("#{blast_dir}/formatdb -p T -i #{non_ambig_dir}/endo_all_nonAmbig.fasta -n #{non_ambig_dir}/endoBlastdb")
end


#Combine all fasta files into 1
def construct_nonAmbig_blastdb(blast_dir)
  non_ambig_dir = "endo_nonAmbig"
  proteins=Dir.glob("#{non_ambig_dir}/*/seq_aa.fasta")
  fastas_to_combine = ""
  proteins.each do |f|
    fastas_to_combine = "#{fastas_to_combine} #{f}"
  end
  system ("cat #{fastas_to_combine} > #{non_ambig_dir}/endo_all_nonAmbig.fasta\n")
  system ("#{blast_dir}/formatdb -p T -i #{non_ambig_dir}/endo_all_nonAmbig.fasta -n #{non_ambig_dir}/endoBlastdb")
end

#blast all ambig cases against nonAmbig database
def blast_each_ambig(basepath)
  ambig_dir = "endo_ambig"
  Dir.chdir(ambig_dir)
  proteins=Dir.glob("*")
  proteins.each do |f|
    Dir.chdir(f)
    if(!File.exists?("seq_aa.fasta.3.psiblast"))
        print("#{basepath}/run-psiblast.pl seq_aa.fasta --db ../../endo_nonAmbig/endoBlastdb --n_rounds 3")
        system("#{basepath}/run-psiblast.pl seq_aa.fasta --db ../../endo_nonAmbig/endoBlastdb --n_rounds 3")
    end
    Dir.chdir("..")
  end
  Dir.chdir("..")
end

def construct_pssm_each_nonAmbig()
  non_ambig_dir = "endo_nonAmbig"
  upSeq_fastas=Dir.glob("#{non_ambig_dir}/*/upSeq_target_dna.fasta")
  downSeq_fastas=Dir.glob("#{non_ambig_dir}/*/downSeq_target_dna.fasta")
  upSeq_fastas.each do |f|
    construct_pssm(f)
  end
  downSeq_fastas.each do |f|
    construct_pssm(f)
  end
end

#used to convert dna seq to fasta file
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


def identify_target_each_ambig(basepath)
  ambig_dir = "endo_ambig"
  numb_blast_hits = 3
  top_hits = Dir.glob("#{ambig_dir}/*/seq_aa.fasta.3.psiblast")
  e_value = BigDecimal.new("9E999") #Very high e_value so all items are read
  top_hits.each do |f|
    blast_results_hash = read_psiblast(f,e_value)
    upSeq_target_dna = "#{File.dirname(f)}/upSeq_target_dna.fasta"
    downSeq_target_dna = "#{File.dirname(f)}/downSeq_target_dna.fasta"
    #get all keys. Put in array. Sort array. If top key print out A. If top 3 print out Y
    blast_e_vals = Array.new
    blast_results_hash.each_value do |blast_hash_v|
        blast_e_val = blast_hash_v.e_val
        blast_e_vals.push(blast_e_val)
    end
    blast_e_vals.sort!
    top_eval = blast_e_vals[0]
    top_3_eval = blast_e_vals[blast_e_vals.size-1]
    if(blast_e_vals.size > 3)
        top_3_eval = blast_e_vals[2]
    end
    blast_results_hash.each_value do |blast_hash_v|
      blast_name =  blast_hash_v.name
      blast_e_val =  blast_hash_v.e_val
      ct = 1
      while((ct <= 3) && (ct <= blast_results_hash.size))
          ct+=1
        if(blast_e_val <= top_3_eval)
            upSeq_dna_fasta = "#{File.dirname(f)}/upSeq_dna.fasta"
            downSeq_dna_fasta = "#{File.dirname(f)}/downSeq_dna.fasta"
            upSeq_template_pssm = "#{blast_name}/upSeq_target_dna.fasta.pssm"
            downSeq_template_pssm = "#{blast_name}/downSeq_target_dna.fasta.pssm"
            upSeq_out = "#{File.dirname(f)}/upSeq_#{File.basename(blast_name)}_#{blast_e_val}.out"
            downSeq_out = "#{File.dirname(f)}/downSeq_#{File.basename(blast_name)}_#{blast_e_val}.out"
        #print("#{basepath}/pssm -s #{upSeq_dna_fasta} -p #{upSeq_template_pssm} -mute > #{upSeq_out}\n")
            system("#{basepath}/pssm -s #{upSeq_dna_fasta} -p #{upSeq_template_pssm} -mute > #{upSeq_out}")
        #print("#{basepath}/pssm -s #{downSeq_dna_fasta} -p #{downSeq_template_pssm} -mute > #{downSeq_out}\n")
            system("#{basepath}/pssm -s #{downSeq_dna_fasta} -p #{downSeq_template_pssm} -mute > #{downSeq_out}")
        end
      end
      if(blast_e_val <= top_eval)#best case
        fl_top_e_val = File.open("#{File.dirname(f)}/top_e_val","w")
        fl_top_e_val << "#{top_eval}\n"
        fl_top_e_val.close
        upSeq_out = "#{File.dirname(f)}/upSeq_#{File.basename(blast_name)}_#{blast_e_val}.out"
        downSeq_out = "#{File.dirname(f)}/downSeq_#{File.basename(blast_name)}_#{blast_e_val}.out"
        upSeq_target_dna_fasta = "upSeq_target_dna.fasta"
        downSeq_target_dna_fasta = "downSeq_target_dna.fasta"
        extract_pred_target_DNA(downSeq_out,"#{File.dirname(f)}/#{downSeq_target_dna_fasta}")
        extract_pred_target_DNA(upSeq_out,"#{File.dirname(f)}/#{upSeq_target_dna_fasta}")
      end
    end
  end
end

#extracts files from predicted pssm file output and puts it in the output
def extract_pred_target_DNA(pred_file_name,targetDna_file_name)
  pred_file = File.new(pred_file_name,'r')
  last_seq = ""
  proteinName = ""
  score = ""
  while(line = pred_file.gets())
    position = line.split(" ")
    if(position.size()>2)
      if(position[2][0..3] ==">ref")
        proteinName = position[2]
        last_seq = position[1]
        score = position[0]
      end
    end
  end
  targetDna_file = File.new(targetDna_file_name,"w")
  targetDna_file << "#{proteinName}\n"
  targetDna_file << "#{score}\n"
  targetDna_file << "#{last_seq}\n"
  targetDna_file.close()
end

#generate two trees, one upsite one down
def generate_cutSite_labeled_trees(basepath)
  rename_all_fasta()
  concat_all_fasta()
  make_phylogenic_trees(basepath)
  #convert_newick_to_nexus("downSeq.tree","downSeq.nex")
  convert_newick_to_nexus("all_upSeq_aa.tree","all_upSeq_aa.nex")
  convert_newick_to_nexus("all_downSeq_aa.tree","all_downSeq_aa.nex")
end

#convert tree from newick to nexus adding the correct coloring
def convert_newick_to_nexus(newick_name,nexus_name)
  newick_file = File.open(newick_name,"r")
  nexus_file = File.open(nexus_name,"w")
  items = Array.new
  while(line = newick_file.gets())
    position = line.split("|")
    if (position.size() > 1)
      items.push(position[0])
    end
  end
  nexus_file << "#NEXUS\n"
  nexus_file << "begin taxa;\n"
  nexus_file << "\tdimensions ntax=#{items.size()};\n"
  nexus_file << "\ttaxlabels\n"
  items.each do |f|
    nexus_file << "\t'#{f}'\n"
  end
  nexus_file << ";\n"
  nexus_file << "end;\n"
  nexus_file << "begin trees;\n"
  nexus_file << "\ttree tree_1 = [&R]\n"
  newick_file.rewind
  while(line = newick_file.gets())
    position = line.split("|")
    if (position.size() > 1)
      nexus_file << "'#{position[0]}'"
      if(position[1] == "known")
        nexus_file << "[&nodeColor=\"255,0,0\"]"
      else if (position[1] == "pred")
             nexus_file << "[&nodeColor=\"0,255,0\"]"

           else
             #nexus_file << "[&nodeColor=\"0,0,255\"]"
           end
      end
      nexus_file << position[2]
    else
      nexus_file << line
    end
  end
  nexus_file << "end;"
end


def make_phylogenic_trees(basepath)
  system("#{basepath}/muscle -in all_downSeq_aa.fasta -tree1 all_downSeq_aa.tree")
  system("#{basepath}/muscle -in all_upSeq_aa.fasta -tree1 all_upSeq_aa.tree")
end

#concatenate all the fasta together
def concat_all_fasta()
  file_type = "aa"
  all_downSeq_files = Array.new
  all_downSeq_files.concat(Dir.glob("endo_known/*/downSeq_#{file_type}.fasta"))
  all_downSeq_files.concat(Dir.glob("endo_ambig/*/downSeq_#{file_type}.fasta"))
  all_downSeq_files.concat(Dir.glob("endo_nonAmbig/*/downSeq_#{file_type}.fasta"))
  all_upSeq_files = Array.new
  all_upSeq_files.concat(Dir.glob("endo_known/*/upSeq_#{file_type}.fasta"))
  all_upSeq_files.concat(Dir.glob("endo_ambig/*/upSeq_#{file_type}.fasta"))
  all_upSeq_files.concat(Dir.glob("endo_nonAmbig/*/upSeq_#{file_type}.fasta"))
  all_downSeq_file_string = ""
  all_downSeq_files.each do |f|
    all_downSeq_file_string.concat("#{f} ")
  end
  all_upSeq_file_string = ""
  all_upSeq_files.each do |f|
    all_upSeq_file_string.concat("#{f} ")
  end
  system("cat #{all_downSeq_file_string} > all_downSeq_aa.fasta\n\n")
  system("cat #{all_upSeq_file_string} > all_upSeq_aa.fasta\n\n")
end


#renames all the fasta so as to include the cutsite in the name
def rename_all_fasta()
  all_endos = Dir.glob("endo_known/*")
  all_endos.concat(Dir.glob("endo_ambig/*"))
  all_endos.concat(Dir.glob("endo_nonAmbig/*"))
  #fix known endos
  #known_endos.each do |f|
  #  system("cp #{f}/#{File.basename(f)}.fasta #{f}/seq_aa.fasta\n")
  #end
  all_endos.each do |f|
    if(!File.file?(f))
      downSeq_cutSite = get_fasta_seq("#{f}/downSeq_target_dna.fasta").downcase
      upSeq_cutSite = get_fasta_seq("#{f}/upSeq_target_dna.fasta").downcase
      #upSeq_orig_fileName = "#{f}/seq_aa.fasta"
      #downSeq_orig_fileName = "#{f}/seq_aa.fasta"
      upSeq_orig_fileName = "#{f}/upSeq_target_dna.fasta"
      downSeq_orig_fileName = "#{f}/downSeq_target_dna.fasta"
      downSeq_fileName = "#{f}/downSeq_aa.fasta"
      upSeq_fileName = "#{f}/upSeq_aa.fasta"
      endo_type = ""
      if(f[5..9] == "known")
        endo_type = "known"
      end
      if(f[5..9] == "ambig")
        endo_type = "ambig"
      end
      if(f[5..9] == "nonAm")
        endo_type = "pred"
      end
      downSeq_internalName = ""
      upSeq_internalName = ""
      if endo_type == "known"
        downSeq_internalName = "#{File.basename(f)}.#{downSeq_cutSite}|#{endo_type}|"
        upSeq_internalName = "#{File.basename(f)}.#{upSeq_cutSite}|#{endo_type}|"
      else
        fullName = File.basename(f)
        partName = fullName.split(".")
        shortname = partName[partName.size()-2]
        downSeq_internalName = "#{shortname}.#{downSeq_cutSite}|#{endo_type}|"
        upSeq_internalName = "#{shortname}.#{upSeq_cutSite}|#{endo_type}|"
      end
      rename_fasta(downSeq_orig_fileName,downSeq_fileName,downSeq_internalName)
      rename_fasta(upSeq_orig_fileName,upSeq_fileName,upSeq_internalName)
    end
  end
end
#|ambig|
#|pred|


#change name of fasta
def rename_fasta(orig_fileName,new_fileName,internalName)
  seq = get_fasta_seq(orig_fileName)
  fasta_file = File.open(new_fileName,"w")
  fasta_file << ">#{internalName}\n"
  fasta_file << "#{seq}\n"
  fasta_file.close()
end

#get content from fasta(my own format?)
def get_fasta_seq(fasta_file_name)
  fasta_file = File.open(fasta_file_name,"r")
  line = fasta_file.gets()
  fasta_seq = ""
  while(line = fasta_file.gets())
    fasta_seq.concat(line.chomp)
  end
  return fasta_seq
end
