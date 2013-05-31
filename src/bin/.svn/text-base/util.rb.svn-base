#!/usr/bin/ruby

#read all config data in
def readConfig(config_file)
  configHash = Hash.new
  File.open(config_file) do |f|
    f.each do |line|
      position = line.split(" ")
      configHash[position[0]]=position[1]
    end
  end
  return configHash
end
