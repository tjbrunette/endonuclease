#!/usr/bin/ruby
require "optparse"
require "set"
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
  opts.on('','--in_fasta FILE(s)',Array,"location of input fastas(comma seperated):REQ") do |f|
    options[:in_fasta] = f
  end
  opts.on('','--in_fasta_dir FILE', String, "location of fasta files") do |f|
    options[:in_fasta_dir] = f
  end
  opts.on( '-h','--help', 'Display this screen') do
    puts opts
    exit
  end
end
optparse.parse!
config_file = options[:config_file]
in_fastas = options[:in_fasta]
if(in_fastas.size() == 0) #read fastas from dir
  in_fasta_dir = options[:in_fasta_dir]
  in_fastas = Dir.glob("#{options[:in_fasta_dir]}/*.fasta")
end
basepath = options[:basepath]
require "#{basepath}/util"
require "#{basepath}/psiblast_util"
require "#{basepath}/gbff_util"
require "#{basepath}/process_endo_results"
configHash = readConfig(config_file)
blastp_db = configHash["blastp_db"]
e_value = configHash["e_value"]
n_rounds_blast = configHash["n_rounds_blast"]
gbff_dir = configHash["gbff_dir"]
seq_dir = configHash["seq_dir"]
flanking_res = configHash["flanking_res"]
blast_dir = configHash["blast_dir"]
psipred_dir = configHash["psipred_dir"]
all_blast_results = Array.new
#/Step 1: Run psiblast on all input fastas
#run_psiblast(in_fastas,blastp_db,e_value,basepath,n_rounds_blast)
#/Step 2: Gather the results from the blast runs
in_fastas.each do |fasta|
  tmp_blast_results = read_psiblast("#{File.basename(fasta)}.#{n_rounds_blast}.psiblast")
  all_blast_results = all_blast_results.concat(tmp_blast_results)
  all_blast_results.each do |f|
    print("#{f.name} #{f.parent} #{f.e_val}\n")
  end
end
#/Step 3: gather all endonucleases
#gather_target_sites(all_blast_results,gbff_dir,seq_dir,flanking_res)
#/Step 4: make prediction
#match_ambig_2_closestNonAmbig(blast_dir,basepath)
#/Step 5: Make cutsite labeled tree.
#generate_cutSite_labeled_trees(basepath)
