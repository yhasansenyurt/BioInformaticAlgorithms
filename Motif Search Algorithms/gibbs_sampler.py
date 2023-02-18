import random
import timeit

#Calculating each motif's score
def find_score(motifs,k,DNA_lines):
    
    score_list=([0]*k)
    
    for y in range(0,k):
        temp_score = ([0] * 4)#Checking A, T, C and G counts in each column.
        for i in range(0,len(DNA_lines)):#

            if (motifs[i][y]=="A"):
                temp_score[0]=temp_score[0]+1
            if (motifs[i][y]=="T"):
                temp_score[1]=temp_score[1]+1
            if (motifs[i][y]=="C"):
                temp_score[2]=temp_score[2]+1
            if (motifs[i][y]=="G"):
                temp_score[3]=temp_score[3]+1
        
        #Initial value is A's count, if there is more counted value than A, then most_repeated will be it.
        most_repeated=temp_score[0]
        
        for p in temp_score: 
            if p>most_repeated:
                most_repeated=p
        score_list[y]=len(DNA_lines)-most_repeated
    
    return sum(score_list)

#selecting 10 random motifs for initial motif set.
def select_random_motifs(DNA_lines,k):
    rows, cols = (len(DNA_lines), k-1)
    motifs  = [[0] * cols] * rows #2D array is created for each line as K-mer number.
    i=0
    
    for lines in DNA_lines:#FROM EACH LINE IT RANDOM TAKES A STRING THE LENGTH OF K-MER AND WRITES INTO THE MOTIFS SEARCH
                          
        random_position=random.randint(0,len(lines)-k)
        k_mer=lines[random_position:random_position+k]
        motifs[i]=k_mer
        i += 1

    return motifs


def create_profil(motifs,k):
    profil = {"A": [1] * k, "C": [1] * k, "T": [1] * k, "G": [1] * k}# Creating profile dictionary with +1 count because gibbs sampler needs it.
    for z in range(0, k):#Each nucleotid count on profile will be increased according to the search.
        
        for y in motifs:
            if y[z] == "A":
                profil["A"][z] = profil["A"][z] + 1
            if y[z] == "T":
                profil["T"][z] = profil["T"][z] + 1
            if y[z] == "C":
                profil["C"][z] = profil["C"][z] + 1
            if y[z] == "G":
                profil["G"][z] = profil["G"][z] + 1
    
    return profil
    
def select_new_motif(profile,removed_DNA_line,k):

    best_matched_motif = ""
    motif_list = list()
    motif_score_list = list()

    for y in range(0,len(removed_DNA_line)-k):#Searching removed line of DNA.
        match_for_each_motif = 1              #We look match for each motif.
        motif=removed_DNA_line[y:y+k]         #k-mer
        motif_list.append(motif)


        for motif_index in range(0,k):        #Calculating scores of each K-MER

            if motif[motif_index]=="A":
                match_for_each_motif=profile["A"][motif_index]*match_for_each_motif
                
            elif motif[motif_index]=="T":
                match_for_each_motif=profile["T"][motif_index]*match_for_each_motif
                                
            elif motif[motif_index]=="C":
                match_for_each_motif=profile["C"][motif_index]*match_for_each_motif
                
            elif motif[motif_index]=="G":
                match_for_each_motif=profile["G"][motif_index]*match_for_each_motif
                
        motif_score_list.append(match_for_each_motif)   #append each k-mers score to select in future with biased die.

    
    #BIASED DIE PART
    motif_score_list = [x / sum(motif_score_list) for x in motif_score_list] #Biased values.
    index = random.choices(range(0, len(motif_score_list)), weights=(motif_score_list))[0] #Throwing die according to bias values.
    best_matched_motif = motif_list[index] #selecting motif according to biased die result.

    return best_matched_motif

def consensus(motifs,k,DNA_lines):
    consensus= ""
    
    for y in range(0,k):

        temp_score = ([0] * 4)#A T C G will be counted.
        for i in range(0,len(DNA_lines)):#

            if (motifs[i][y]=="A"):
                temp_score[0]=temp_score[0]+1
            if (motifs[i][y]=="T"):
                temp_score[1]=temp_score[1]+1
            if (motifs[i][y]=="C"):
                temp_score[2]=temp_score[2]+1
            if (motifs[i][y]=="G"):
                temp_score[3]=temp_score[3]+1

        most_repeated=temp_score[0]#most repeated value is A for inital value.

        index = 0
        for p in range (0,4):
            if temp_score[p] > most_repeated:
                most_repeated = temp_score[p]
                index = p
        #Consensus is being creating now.
        if index == 0:
            consensus += 'A'
        elif index == 1:
            consensus += 'T'
        elif index == 2:
            consensus += 'C'
        elif index == 3:
            consensus += 'G'
        
    return consensus

def gibbs_sampler(DNA,k):

    DNA_lines = DNA.readlines()
    motifs = select_random_motifs(DNA_lines,k) #random 10 k-mer for the beginning of algorithm.
    bestMotifs = motifs.copy()

    iteration_to_stop = 0   #if it is greater than 100, it means in last 100 iterations, score did not improved. So algorithm will stop.
    iteration_number = 0    #Checking score in each 50 iterations to decide whether algorithm will stop or not.

    while True:
        
        removed_kmer_index = random.randint(0,9)
        motifs.pop(removed_kmer_index)                                  #9 motif.
        removed_DNA_line = DNA_lines[removed_kmer_index]                #Removed DNA line.
        
        profile = create_profil(motifs,k)                               #+1 li profil
        selected_motif = select_new_motif(profile,removed_DNA_line,k)   #new motif is selected.

        motifs.insert(removed_kmer_index,selected_motif)                #10 motif now.
        
        #Calculating scores of motifs to compare them to select best one.
        motif_score = find_score(motifs,k,DNA_lines)
        bestmotif_score = find_score(bestMotifs,k,DNA_lines)
        
        iteration_number +=1

        #Comparing scores.
        if motif_score < bestmotif_score:
            bestMotifs = motifs.copy()
            iteration_to_stop = 0
        else:
            iteration_to_stop +=1
        #Stop or continue.
        if iteration_to_stop >= 100 and iteration_number % 50 == 0:
            break
            
    print("BEST MOTIFS: ",bestMotifs)
    print("SCORE: ",find_score(bestMotifs,k,DNA_lines))
    print("CONSENSUS: ",consensus(bestMotifs,k,DNA_lines))
    print("Iteration number: ",iteration_number)
    return bestMotifs

def main():

    DNA = open("input.txt", "r")
    k = 10

    start = timeit.default_timer()
    gibbs_sampler(DNA,k)
    stop = timeit.default_timer()

    print('Time: ', stop - start)


main()