#!/usr/bin/ruby
require "fileutils"
require "set"
require "bigdecimal"
#require "ftools"
#run psiblast
def run_psiblast(in_fastas,blastp_db,e_value,basepath,n_rounds,blast_ver)
  in_fastas.each do |fasta|
    system("cp #{fasta} .")
    system("#{basepath}/run-psiblast2.pl #{File.basename(fasta)} --db #{blastp_db} --e_value #{e_value} --h_value #{e_value} --n_rounds #{n_rounds} --blastpgp #{blast_ver}\n")
  end
end

#psiblast_result class
class Psiblast_result
  attr_accessor :seq, :parent, :e_val, :name, :output_to_db
  def initialize(name,seq,e_val,parent)
    @seq = seq
    @parent = parent
    @e_val = e_val
    @name = name
    @output_to_db = false
  end
end

#Read and parse a psiblast file into an array Psiblast_result
#it should be noted that the e-values come only from the final round.
def read_psiblast(in_psiblast,max_e_value)
  blastEntries = Hash.new
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
      if((subposition[0] == ">ref")||(subposition[0] == ">gb")||(subposition[0] == ">sb")||(subposition[0] == ">dbj")||(subposition[0] == ">sp")||(subposition[0] == ">emb")||(subposition[0] == ">pdb")||(subposition[0] == ">prf") )
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
          e_val_tmp = position_2nd[7][0..-2]
          #e_val_tmp = (position_2nd[7][0..position_2nd[7].size()-2])
          if(e_val_tmp.to_s[0..0] == "e")
              e_val_tmp = "1#{e_val_tmp}"
          end
            e_val = BigDecimal.new(e_val_tmp)
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
            if(position_2nd[0] != nil)
               if(position_2nd[0].include?("Sbjct"))
                seq = seq + position_2nd[2]
               end
            end
          end
        end
        if((e_val <= max_e_value) && (seq.size >1))
            blast_item = Psiblast_result.new(name,seq,e_val,parent)
            blast_item_array = blast_item.clone
            blastEntries[blast_item.name] = blast_item_array
        end
      end
    end
  end
  return blastEntries
end

def read_blast_hash(blast_hash_fn)
    fl = File.open(blast_hash_fn)
    blastEntries = Hash.new
    while(line = fl.gets)
        pos = line.split(" ")
        name = pos[0]
        seq = pos[4]
        parent =pos[3]
        e_val = pos[2]
        if(seq != nil)
            adjusted_seq = seq.gsub("-","")
            if(adjusted_seq.size > 1)
                blast_item =Psiblast_result.new(name,adjusted_seq,e_val,parent)
                blast_item_array = blast_item.clone
                blastEntries[blast_item.name] = blast_item_array
            end
        end
    end
    return(blastEntries)
end
