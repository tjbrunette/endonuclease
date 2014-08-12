#This is to be run after you have pulled out sheets,
#clustered with muscle, and generated a distance_matrix using the following command
#~krypton/bin/nw_distance -n -m m beta_sheet.tree  > distance_tree
distance_matrix = File.open("distance_tree","r")
line = distance_matrix.gets
names = line.split(" ")
non_ambig_hash = Hash.new
for ii in 0..names.size-1 do
  if(names[ii].split("/")[0] == "endo_nonAmbig")
    non_ambig_hash[ii] = names[ii]
  end
end
while(line = distance_matrix.gets)
  score_array = line.split(" ")
  if(score_array[0].split("/")[0] != "endo_nonAmbig") #No need to match the non_ambiguous
    distance_hash = Hash.new
    non_ambig_hash.each do |key,value|
      distance_hash[score_array[1+key]] = value
    end
    lowest_dist = distance_hash.keys.sort[0]
    closest_ambig = distance_hash[lowest_dist]
    system("cp #{closest_ambig}/upSeq_target_dna.fasta #{score_array[0]}/closest_upSeq_target_dna.fasta\n")
    system("cp #{closest_ambig}/downSeq_target_dna.fasta #{score_array[0]}/closest_downSeq_target_dna.fasta\n")
  end
end
