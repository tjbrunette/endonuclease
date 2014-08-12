#!/usr/bin/ruby

#gathers target site data
def gather_target_sites(all_blast_results,gbff_dir,seq_dir,flanking_res)
  #step1 process blast_results into a hash function which includes the seq_id and an array of blast_items
  #blast_result_hash = blast_results_to_hash(all_blast_results)
  #step2 gather endo from gbff file
  process_all_db(all_blast_results,gbff_dir,seq_dir,flanking_res)
end

#create a hash function which hashes the seq_id into an array of blast_items
def blast_results_to_hash(all_blast_results)
  blast_result_hash = Hash.new
  all_blast_results.each do |f|
    if(blast_result_hash.has_key?(f.name))
      blast_result_hash[f.name].push(f)
    else
      tmp_array = f.clone
      blast_result_hash[f.name] = tmp_array
    end
  end
  return(blast_result_hash)
end

#endonuclease_from gbff
class Endo_info
  attr_accessor :name,:upseq_target,:downseq_target,:upseq,:downseq,:dnaStartRes,:dnaStopRes,:aaSeq,:ambig,:organismName
  def initialize(name,upseq_target,downseq_target,upseq,downseq,dnaStartRes,dnaStopRes,aaSeq,ambig,organismName)
    @name = name
    @upseq_target = upseq_target
    @downseq_target = downseq_target
    @upseq =upseq
    @downseq= downseq
    @dnaStartRes = dnaStartRes
    @dnaStopRes = dnaStopRes
    @aaSeq = aaSeq
    @ambig = ambig
    @organismName = organismName #currently Fred
  end
end

#gene_from gbff
class Gene_info
  attr_accessor :name,:startPos,:stopPos,:aaSeq,:known_endo
  def initialize(name,startPos,stopPos,aaSeq,known_endo)
    @name = name
    @startPos = startPos
    @stopPos = stopPos
    @aaSeq = aaSeq
    @known_endo = known_endo
  end
end
#gene_from gbff
class Exon_info
  attr_accessor :startPos,:stopPos
  def initialize(startPos,stopPos)
    @startPos = startPos
    @stopPos = stopPos
  end
end

def process_single_db(blast_result_hash,flanking_res,gbff_file)
  fna_hash = Hash.new
  fna_file = "#{gbff_file[0..-6]}.fna"
  if(File.exists?(fna_file))
    fna_hash = process_fna(fna_file)
  end
  print("processing gbff #{gbff_file}\n")
  process_gbff(blast_result_hash,gbff_file,flanking_res,fna_hash) #also process seq files
  print("completed processing #{gbff_file}\n")
end

#process all gbff files
def process_all_db(blast_result_hash,gbff_dir,seq_dir,flanking_res)
  gbff_files = Array.new
  numbThreads = 3
  gbff_files =Dir.glob("#{seq_dir}/*.seq")
  gbff_files.concat(Dir.glob("#{gbff_dir}/*.gbff"))
  gbff_files.each do |gbff_file|
    process_single_db(blast_result_hash,flanking_res,gbff_file)
    #print("#{Thread.list.size}\n")
    # if(Thread.list.size < numbThreads)
    #  a = Thread.new{process_single_db(blast_result_hash,flanking_res,gbff_file)}
     # a.run
    #else
    #  Thread.list.each do |thr|
    #    thr.join
    #  end
    #end
  end
  #output all the unseen items in the database to the insufficient DNA database
    blast_result_hash.each do |key,value|
        if(value.output_to_db == false)
            output_blast_hit(value)
        end
    end
end


#process fna file by extracting dna
def process_fna(fna_file)
  fna_hash = Hash.new
  input = File.new(fna_file,'r')
  first_loop = true
  dna_seq = ""
  id = ""
  count = 0
  while(line = input.gets)
    if(line[0..0] == '>')
      if(first_loop == false)
        fna_hash[id] = dna_seq
      end
      position=line.split("|")
      id = position[3][0..-3]
      dna_seq = ""
      first_loop = false
    else
      count = count+1
      dna_seq.concat(line[0..-2])
    end
  end
  fna_hash[id] = dna_seq
  return(fna_hash)
end

#gather all endonucleases from file
def process_gbff(blast_result_hash,gbff_db,flanking_res,fna_hash)
  gbff_file = File.new(gbff_db,'r')
  while(line = gbff_file.gets)
    position = line.split()
    if(position[0] == "LOCUS")
      endo_proteins = process_organism_gbff(gbff_file,blast_result_hash,flanking_res,fna_hash)
      output_endo_proteins(endo_proteins,flanking_res)
    end
  end
end


#process organism
def process_organism_gbff(gbff_file,blast_result_hash,flanking_res,fna_hash)
  endo_proteins = Array.new
  organism_data = Array.new
  organism_data_type = Array.new
  dna_seq = ""
  accession = ""
  organism_done = false
  while(organism_done == false)
    if(gbff_file.eof? == true)
      print("partial organism encountered at end of file.")
      organism_done = true
    else 
        line = gbff_file.gets
        position = line.split()
        if(position[0] == "ORIGIN")
            dna_seq = process_ORIGIN(gbff_file)
            organism_done = true
        end
        if(position[0] == "CONTIG")
            if fna_hash.has_key?(accession)
                dna_seq = fna_hash[accession]
            else
                dna_seq = ""
            end
        organism_done = true
        end
        if(position[0] == "ACCESSION")
            accession = process_ACCESSION(line)
        end
        if(position[0] == "exon")
            if(position.size() > 1)
                organism_data.push(process_exon(line))
                organism_data_type.push("exon")
            end
        end
        if(position[0] == "CDS" && position.size > 1 && (line[0..7] != "        "))
            #lineNumb = gbff_file.lineno
            gene = process_CDS(line,gbff_file)
            if(gene != nil)
                organism_data.push(gene)
                organism_data_type.push("gene")
            else
                #gbff_file.rewind
                #(lineNumb).times{gbff_file.gets}
                line = gbff_file.gets
            end
        end
    end
  end
  for i in 0..(organism_data.size()-1)
    if(organism_data_type[i] == "gene")
      if(blast_result_hash.key?(organism_data[i].name)||(organism_data[i].known_endo))
          #tracks which items from the blast DB have been found.
          if(blast_result_hash.key?(organism_data[i].name))
              blast_result_hash[organism_data[i].name].output_to_db = true
          end
        endo_inserted = false
        if((i-1)>=0 && (i<organism_data_type.size()-1))
          if(organism_data_type[i-1] == "exon" && organism_data_type[i+1] == "exon")
            #non ambiguous case
            endo_proteins.push(process_endo(organism_data[i-1],organism_data[i+1],organism_data[i],dna_seq,accession,flanking_res,false))
            endo_inserted = true
          end
        end
        if(endo_inserted == false)
          #ambiguous case
          endo_proteins.push(process_endo("","",organism_data[i],dna_seq,accession,flanking_res,true))
        end
      end
    end
  end
  return(endo_proteins)
end

#get accession number
def process_ACCESSION(line)
  position=line.split(" ")
  return(position[1])
end

#process dna_seq
def process_ORIGIN(gbff_file)
  dna_seq = ""
  line = gbff_file.gets
  position = line.split(" ")
  while(position[0] != "//")
    for i in (1..position.size()-1)
      dna_seq = dna_seq.concat(position[i])
    end
    line = gbff_file.gets
    position = line.split(" ")
  end
  return(dna_seq)
end

#process exon
def process_exon(line)
  position = line.split(" ")
  exonString = position[1]
  if(exonString[0..9] == "complement")
    exonString = exonString[11..exonString.size()-2]
  end
  subPosition = exonString.split("..")
  startPos = subPosition[0].to_i
  stopPos = subPosition[1].to_i
  exon = Exon_info.new(startPos,stopPos)
end

def process_CDS(line,gbff_file)
  position = line.split(" ")
  genePosString = position[1]
  startPos = -1
  endPos = -1
  known_endo = false
  if(genePosString[0..3] == "join")
    genePosString = genePosString[5..-2]
  end
  if(genePosString[0..9] == "complement")
    genePosString = genePosString[11..-2]
  end
  if(genePosString[0..3] == "join")
    genePosString = genePosString[5..-2]
  end
  if(genePosString.include?(":"))
    genePosString = genePosString.split(":")[1]
  end
  #process the rest of the string
  subPosition = genePosString.split("..")
  if(subPosition.size == 1)
    startPos = subPosition[0].to_i
    stopPos = subPosition[0].to_i
  else
    if(genePosString.include?(","))
      tmp = subPosition[1].split(",")
      subPosition[1] = tmp[0]
    end
    if(subPosition[0][0..0] == "<")
      startPos = subPosition[0][1..-1].to_i
    else
      startPos = subPosition[0].to_i
    end
    if(subPosition[1][0..0] == ">")
      stopPos = subPosition[1][1..-1].to_i
    else
      stopPos = subPosition[1].to_i
    end
  end
  #process each part of the gbff file
  if(gbff_file.eof? == true)
    print "dieing here***A\n"
    Process.exit()
  end
  line = gbff_file.gets
  len = line.size
  position = line.split("=")
  bad_CDS_position = line.split(" ")
 # while((position[0].lstrip != "/translation") && (bad_CDS_position[0] != "gene") && (bad_CDS_position[0] != "misc_feature") && (bad_CDS_position[0] != "mRNA")&& (bad_CDS_position[0] != "exon") && (bad_CDS_position[0] != "FEATURES") && (bad_CDS_position[0] != "source") && (bad_CDS_position[0] != "/db_xref=\"PSEUDO:CBN59414.1\""))
  proteinid = ""
  aaSeq =""
  while((line[0..7] == "        ") && (position[0].lstrip != "/translation"))
    if(position[0].lstrip == "/protein_id")
      proteinid = position[1][1..position[1].size-3]
    end
    if((position[0].lstrip == "/product")||(position[0].lstrip == "/note"))
        if(position[1].include?("LAGLIDAGD")||position[1].include?("LAGLIDADG")) 
      #if(position[1][1..9] == "LAGLIDADG"||position[1][1..9] == "LAGLIDAGD")
        known_endo = true
      end
    end
    if(gbff_file.eof? == true)
      print "dieing here***B\n"
      Process.exit()
    end
    line = gbff_file.gets
    len += line.size
    position = line.split("=")
    bad_CDS_position = line.split(" ")
  end
 # if((bad_CDS_position[0] != "gene")&&(bad_CDS_position[0] != "misc_feature")&&(bad_CDS_position[0] != "mRNA")&&(bad_CDS_position[0] != "exon") && (bad_CDS_position[0] != "FEATURES") && (bad_CDS_position[0] != "source") && (bad_CDS_position[0] != "/db_xref=\"PSEUDO:CBN59414.1\""))
    #process aa seq
  if(position[0].lstrip == "/translation")
    lineWoSpace = position[1].lstrip[1..-1]#remove first "
    while(lineWoSpace[-2..-2] != "\"")
      aaSeq = aaSeq.concat(lineWoSpace[0..-2])
      line = gbff_file.gets
      lineWoSpace =line.lstrip
    end
    aaSeq = aaSeq.concat(lineWoSpace[0..-3])
  end
  gene = nil
  if((!proteinid.empty?)&&(!aaSeq.empty?))
    gene = Gene_info.new(proteinid,startPos,stopPos,aaSeq,known_endo)
  else
    gbff_file.seek(-len,IO::SEEK_CUR)
  end
  return(gene)
end

#process non-ambiguous endonucleases
def process_endo(preExon,postExon,gene,dna_seq,accession,flanking_res,ambig)
  exonSize = 15
  #**DNA is stored 0 indexed but the info in the gbff is 1 indexed
  #first get the predicted target.
  upSeq_target = ""
  downSeq_target = ""
  if(ambig == false)
    upSeq_start =preExon.stopPos.to_i-exonSize-1
    upSeq_end = preExon.stopPos.to_i-1
    downSeq_start = postExon.startPos.to_i-1
    downSeq_end = postExon.startPos.to_i+exonSize-1
    upSeq_target = dna_seq[upSeq_start..upSeq_end]
    downSeq_target = dna_seq[downSeq_start..downSeq_end]
  end
  #second get the flanking res
  upSeq_start =gene.startPos.to_i-flanking_res.to_i-1
  upSeq_end = gene.startPos.to_i-1
  downSeq_start = gene.stopPos.to_i-1
  downSeq_end = gene.stopPos.to_i+flanking_res.to_i-1
  if(upSeq_start < 0) #the case where the gene is within flanking_res to the beginning of the chain
    upSeq_start = 0
  end
  if(downSeq_end > (dna_seq.size()-1))#the case where the flanking_res would stretch beyond the end of the gene
    downSeq_end = dna_seq.size()-1
  end
  upSeq = ""
  downSeq = ""
  if((dna_seq.size != 0) && (upSeq_start<dna_seq.size) && (downSeq_start < dna_seq.size))
    upSeq = dna_seq[upSeq_start..upSeq_end]
    downSeq = dna_seq[downSeq_start..downSeq_end]
  end
  endo = Endo_info.new(gene.name,upSeq_target,downSeq_target,upSeq,downSeq,upSeq_start+1,downSeq_end+1,gene.aaSeq,ambig,accession)
  return(endo)
end


#output endonuclease proteins
def output_endo_proteins(endo_proteins,flanking_res)
  if(!File.exists?("endo_ambig"))
    Dir.mkdir("endo_ambig")
  end
  if(!File.exists?("endo_nonAmbig"))
    Dir.mkdir("endo_nonAmbig")
  end
  if(!File.exists?("endo_insufficient_DNA"))
    Dir.mkdir("endo_insufficient_DNA")
  end
  endo_proteins.each do |endo|
    output_endo(endo,flanking_res)
  end
end
#output blast hit
def output_blast_hit(blast_hit)
    endo = Endo_info.new(blast_hit.name,"","","","","0","0",blast_hit.seq.gsub("-",""),true,"blastFound")
    output_endo(endo,"")
end

#output endo
def output_endo(endo,flanking_res)
  all_out = "all_info.txt"
  upSeq_dna = "upSeq_dna.fasta"
  downSeq_dna = "downSeq_dna.fasta"
  upSeq_target_dna = "upSeq_target_dna.fasta"
  downSeq_target_dna = "downSeq_target_dna.fasta"
  seq_aa = "seq_aa.fasta"
  dir_name = ""
  endo_name = "#{endo.organismName}.#{endo.dnaStartRes}_#{endo.dnaStopRes}.#{endo.name}"
  insufficientDNA = false
  if endo.ambig == true
    dir_name = "endo_ambig/#{endo_name}"
    if((endo.upseq.size()<flanking_res.to_i)||(endo.downseq.size()<flanking_res.to_i)||flanking_res.to_i == 0)
      dir_name = "endo_insufficient_DNA/#{endo_name}"
      insufficientDNA = true
    end
  else
    dir_name = "endo_nonAmbig/#{endo_name}"
  end
  tmp_filename = "#{dir_name}/#{all_out}"
  if(File.exists?(tmp_filename))
  else
    Dir.mkdir(dir_name)
    output = File.new(tmp_filename, "w")
    output << "name:#{endo.name}\n"
    if(endo.ambig == false)
      output << "upseq_target:#{endo.upseq_target}\n"
      output << "downseq_target:#{endo.downseq_target}\n"
    end
    if(!insufficientDNA)
      output << "upseq:#{endo.upseq}\n"
      output << "downseq:#{endo.downseq}\n"
      output << "dnaStartRes:#{endo.dnaStartRes}\n"
      output << "dnaStopRes:#{endo.dnaStopRes}\n"
    end
    output << "aa:#{endo.aaSeq}\n"
    tmp_filename = "#{dir_name}/#{upSeq_dna}"
    output_fasta(tmp_filename,endo_name,endo.upseq)
    tmp_filename = "#{dir_name}/#{downSeq_dna}"
    output_fasta(tmp_filename,endo_name,endo.downseq)
    if(endo.ambig == false)
      tmp_filename = "#{dir_name}/#{upSeq_target_dna}"
      output_fasta(tmp_filename,endo_name,endo.upseq_target)
      tmp_filename = "#{dir_name}/#{downSeq_target_dna}"
      output_fasta(tmp_filename,endo_name,endo.downseq_target)
    end
    tmp_filename = "#{dir_name}/#{seq_aa}"
    output_fasta(tmp_filename,endo_name,endo.aaSeq)
    output.close()
  end
end

#output fasta
def output_fasta(filename,endo_name,seq)
  output = File.new(filename,"w")
  output << ">ref|#{endo_name}\n"
  output << "#{seq}\n"
  output.close()
end

