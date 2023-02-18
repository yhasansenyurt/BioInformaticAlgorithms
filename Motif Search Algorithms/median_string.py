import itertools
import timeit

#Creating all possible k-mers.
def create_k_mers(k):
    nucleotides = ['A','T','C','G']
    all_k_mers_tuple = itertools.product(nucleotides, repeat=k)
    all_k_mers_str = []
    for j in all_k_mers_tuple:
        all_k_mers_str.append(''.join(j))
    
    return all_k_mers_str #returns list of k-mers

#calculating hamming distance between two strings
def hamming_distance(str1,str2):
    count = 0
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            count +=1
    return count

#searching one k_mer in each line of DNA and calculate total distance of it.
def distance(k_mer, DNA):
    distance = 0
    
    for line in range(0,10):
        line_distance = len(k_mer) #max distance
        DNA_line = DNA[line][:-1]  #each line of DNA
        
        #calculating distance in one line shifting right one by one.
        for i in range(len(DNA_line)-len(k_mer)+1):
            if hamming_distance(k_mer,DNA_line[i:i+len(k_mer)]) < line_distance:
                line_distance = hamming_distance(k_mer,DNA_line[i:i+len(k_mer)])
                #print(DNA_line[i:i+len(k_mer)])
        #print(line_distance)        
        distance += line_distance

    return distance

#median string algorithm
def median_string(k,DNA):
    k_mers = create_k_mers(k)
    best_k_mer = k_mers[0] #first k-mer as default for best k-mer
    
    for k_mer in k_mers:
        #searching for the best k-mer by looking distance values.
        if distance(k_mer,DNA) < distance(best_k_mer,DNA):
            best_k_mer = k_mer
            print(best_k_mer)

    return best_k_mer


def main():
    
    DNA = open("input.txt", "r")
    k = 10
    DNA_lines = DNA.readlines()
    
    start = timeit.default_timer()
    median_string(k,DNA_lines)
    stop = timeit.default_timer()

    print('Time: ', stop - start)  
    

main()