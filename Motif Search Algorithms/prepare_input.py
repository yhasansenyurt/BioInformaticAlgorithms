import random

file = open("input.txt", "w")
nucleotides = ['A','T','C','G']
line = ""
ten_mer = "ATCGGCTATT"


for j in range (0,10):
    for i in range(0,500):
        line += nucleotides[random.randint(0,3)]
    
    #Chosing positions of mutations in ten_mer
    mutation_positions = random.sample(range(10), 4)
    ten_mer = "ATCGGCTATT"

    #adding mutations
    for mp in mutation_positions:
        
        #choosing new and different nucleotid for selected position of 10-mer
        new_nucleotid = nucleotides[random.randint(0,3)]
        while new_nucleotid == ten_mer[mp]:
            new_nucleotid = nucleotides[random.randint(0,3)]
        
        #replacing it.
        ten_mer = ten_mer[:mp] + new_nucleotid + ten_mer[mp+1:]
            
    #inserting 10-mer randomly to each line.
    index = random.randint(0,489)
    print(ten_mer,index)
    line = line[:index] + ten_mer + line[index+10:]
    
    #writing a line to a input file
    file.write(line)
    file.write("\n")
    line = ""


