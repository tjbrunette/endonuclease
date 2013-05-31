endos = Dir.glob("endo_nonAmbig/*")
endos += Dir.glob("endo_insufficient_DNA/*")
endos += Dir.glob("endo_ambig/*")
ct = 0
list = Array.new
for n in 0..19 do
  list.push(File.open("list#{n}","w"))
end
endos.each do |endo|
  list[ct%20] << "#{endo}\n"
  ct+=1
end
list.each do |item|
  item.close
end

