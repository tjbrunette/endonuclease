#!/usr/bin/ruby
require "set"
require "ftools"
#run psiblast
def run_psiblast(in_fastas,blastp_db,e_value,basepath,n_rounds)
  in_fastas.each do |fasta|
    system("cp #{fasta} .")
    system("#{basepath}/run-psiblast.pl #{File.basename(fasta)} --db #{blastp_db} --e_value #{e_value} --n_rounds #{n_rounds}\n")
  end
end

#psiblast_result class
class Psiblast_result
  attr_accessor :seq, :parent, :e_val, :name
  def initialize(name,seq,e_val,parent)
    @seq = seq
    @parent = parent
    @e_val = e_val
    @name = name
  end
end

#Read and parse a psiblast file into an array Psiblast_result
#it should be noted that the e-values come only from the final round.
def read_psiblast(in_psiblast)
  blastEntries = Array.new
  parent = ""
  if(!File.exists?(in_psiblast))
    abort("ERROR Blast File #{in_psiblast} does not exist.\n")
  end
  parent = File.basename(in_psiblast)[0..4]
  file = File.new(in_psiblast,'r')
  count = 0
  while(line = file.gets)
    name = ""
    seq = ""
    e_val = 0
    count=count+1
    position = line.split(" ")
    if(position[0] != nil)
      #get name --
      if(position[2] == "round")
        blastEntries.clear()
      end
      subposition = position[0].split("|")
      if((subposition[0] == ">ref")||(subposition[0] == ">gb")||(subposition[0] == ">sb")||(subposition[0] == ">dbj"))
        name = subposition[1]
        #get e-value
        line = file.gets
        position_2nd = line.split(" ")
        while(position_2nd[0] != "Score")
          line=file.gets
          position_2nd = line.split(" ")
        end
        position_2nd = line.split(" ")
        if(position_2nd[0] == "Score")
          e_val = position_2nd[7][0..position_2nd[7].size()-2]
        end
        #get seq -- ends when 2 blank lines in a row are encountered
        blank_count=0
        line=file.gets
        while(blank_count < 2)
          line=file.gets
          if(line.size() == 1)
            blank_count = blank_count+1
          else
            blank_count = 0
            position_2nd = line.split(" ")
            if(position_2nd[0] == "Sbjct:")
              seq = seq + position_2nd[2]
            end
          end
        end
        blast_item = Psiblast_result.new(name,seq,e_val,parent)
        blastEntries.push(blast_item)
      end
    end
  end
  return blastEntries
end



