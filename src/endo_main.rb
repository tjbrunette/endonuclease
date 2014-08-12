#!/usr/bin/ruby
require "fileutils"
require "optparse"
require "set"
require "bigdecimal"
options = {}
optparse = OptionParser.new do |opts|
  options[:basepath] = "/work/brunette/src/endo/bin"
  options[:config_file] = "/work/brunette/src/endo/endo.config"
  options[:in_fasta] = ""
  options[:in_fasta_dir] = ""
  opts.on('','--config_file FILE',String,"location of config file") do |f|
    options[:config_file] = f
  end
  opts.on('','--basepath FILE',String,"location of executables") do |f|
    options[:basepath] = f
  end
  opts.on('','--fasta FILE(s)',Array,"location of input fastas(comma seperated)") do |f|
    options[:fasta] = f
  end
  opts.on('','--fasta_dir FILE', String, "location of fasta files") do |f|
    options[:fasta_dir] = f
  end
   opts.on('','--blast_dir FILE', String, "location of blast files") do |f|
    options[:blast_dir] = f
  end
    opts.on('','--blast_hash FILE', String, "location of blast hash file") do |f|
    options[:blast_hash] = f
  end

  opts.on( '-h','--help', 'Display this screen') do
    puts opts
    exit
  end
end
optparse.parse!
config_file = options[:config_file]
basepath = options[:basepath]
require "#{basepath}/util"
require "#{basepath}/psiblast_util"
require "#{basepath}/gbff_util"
require "#{basepath}/process_endo_results"
configHash = readConfig(config_file)
blastp_db = configHash["blastp_db"]
e_value = BigDecimal.new(configHash["e_value"])
n_rounds_blast = configHash["n_rounds_blast"]
gbff_dir = configHash["gbff_dir"]
seq_dir = configHash["seq_dir"]
flanking_res = configHash["flanking_res"]
blast_ver = configHash["blast_ver"]
blast_dir = configHash["blast_dir"]
psipred_dir = configHash["psipred_dir"]
all_blast_results = Hash.new
#/Step 1: Run psiblast on all input fastas.  Only run if fasta files are defined
if(options[:fasta] != nil || options[:fasta_dir] != nil)
  if(options[:fasta] != nil)
    in_fastas = options[:fasta]
  else
    #read fastas from dir
    in_fasta_dir = options[:fasta_dir]
    in_fastas = Dir.glob("#{options[:fasta_dir]}/*.fasta")
  end
  run_psiblast(in_fastas,blastp_db,e_value,basepath,n_rounds_blast,blast_ver)
end
#/Step 2: Gather the results from the blast runs
if(options[:blast_hash]!= nil)
    all_blast_results = read_blast_hash(options[:blast_hash])
else
    options[:blast_dir]
    blast_files = Dir.glob("#{options[:blast_dir]}/*.psiblast")
    ct = 0
    blast_files.each do |blast_file|
    print("Working on #{ct} #{blast_file}\n")
    ct+=1
    tmp_blast_results = read_psiblast(blast_file,e_value)
    all_blast_results = all_blast_results.merge!(tmp_blast_results)
    end
end
#temporary-----------------------
name_h= Hash.new
ns = File.open("name_seq.txt")
while(line = ns.gets)
    pos = line.split(" ")
    name_h[pos[0]] = pos[1]
end
#--------------Get all fasta & accession names. 
#blast_hash = File.open("blast_hash.txt","w")
subset_blast_results = Hash.new 
all_blast_results.each do |key,blast_hit|
    if(!name_h.key?(blast_hit.name)) #if you haven't done it before
        #---------also temporary-& first end below----------
        print("endo_insufficient_DNA/blastFound.0_0.#{blast_hit.name}\n")
        subset_blast_results[key] = blast_hit       
    end
end
exit(1)
#/Step 3: gather all endonucleases
gather_target_sites(all_blast_results,gbff_dir,seq_dir,flanking_res)
#/Step 4: make prediction
match_ambig_2_closestNonAmbig(blast_dir,basepath)
#/Step 5: Make cutsite labeled tree.
#generate_cutSite_labeled_trees(basepath)
