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


convert_newick_to_nexus("beta_sheet.tree","beta_sheet.nex")
