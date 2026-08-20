

script_location = try
    ARGS[1]
catch
    "local"
end

if script_location == "HPC"
    println(script_location)
    flag = "--project=~/julia_coding/chapter_one/" # HPC flag
elseif script_location == "local"
    println(script_location)
    flag = "--project=~/Documents/julia_coding/chapter_one/" # local machine flag
end

using Arrow, DataFrames

df = DataFrame(entry = 1)

file = string("../../../../mnt/parscratch/users/bi1tma/test.arrow")

Arrow.write(file, df)