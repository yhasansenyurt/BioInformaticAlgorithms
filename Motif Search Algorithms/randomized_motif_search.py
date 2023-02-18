import timeit
import random


def find_score(motifs,k,DNA_lines):#Gönderilen Motifin skorunu hesaplıyor / score'u integer olarak return ediyor
    score_list=([0]*k)#Her bir sütunun scorunu tutuyor
    
    for y in range(0,k):
        temp_score = ([0] * 4)#A T C G nin her bir sütunda kaç kere olduğuna bakmamız için array
        
        for i in range(0,len(DNA_lines)):#

            if (motifs[i][y]=="A"):
                temp_score[0]=temp_score[0]+1
            if (motifs[i][y]=="T"):
                temp_score[1]=temp_score[1]+1
            if (motifs[i][y]=="C"):
                temp_score[2]=temp_score[2]+1
            if (motifs[i][y]=="G"):
                temp_score[3]=temp_score[3]+1
                
        most_repeated=temp_score[0]#BAŞLANGIÇ DEĞERİ OLARAK A'NIN KAÇ KERE TEKRARLANDIĞINI VERİYORUZ ONDAN FAZLA DEĞER VARSA MOST_REPEATED'A ATIYORUZ
        
        for p in temp_score:
            if p>most_repeated:
                most_repeated=p
        score_list[y]=len(DNA_lines)-most_repeated

    return sum(score_list)

def consensus(motifs,k,DNA_lines):
    consensus= ""
    
    for y in range(0,k):

        temp_score = ([0] * 4)#A T C G nin her bir sütunda kaç kere olduğuna bakmamız için array

        for i in range(0,len(DNA_lines)):#
            if (motifs[i][y]=="A"):
                temp_score[0]=temp_score[0]+1
            if (motifs[i][y]=="T"):
                temp_score[1]=temp_score[1]+1
            if (motifs[i][y]=="C"):
                temp_score[2]=temp_score[2]+1
            if (motifs[i][y]=="G"):
                temp_score[3]=temp_score[3]+1
        most_repeated=temp_score[0]#BAŞLANGIÇ DEĞERİ OLARAK A'NIN KAÇ KERE TEKRARLANDIĞINI VERİYORUZ ONDAN FAZLA DEĞER VARSA MOST_REPEATED'A ATIYORUZ

        index = 0
        for p in range (0,4):
            if temp_score[p] > most_repeated:
                most_repeated = temp_score[p]
                index = p
        if index == 0:
            consensus += 'A'
        elif index == 1:
            consensus += 'T'
        elif index == 2:
            consensus += 'C'
        elif index == 3:
            consensus += 'G'
        
    return consensus


def create_profil(motifs,k,DNA_lines):#Verilen motif arrayine göre profil oluşturuyor
    profil = {"A": [0] * k, "C": [0] * k, "T": [0] * k, "G": [0] * k}# HER BİR NÜKLEOTİT İÇİN DİCTİONARY İÇİNDE K-MER UZUNLUĞUNDA ARRAY OLUŞTURUYOR
    
    for z in range(0, k):#MOTİFTEKİ NÜKLEOTİTLERE GÖRE DİCT İÇİNDEKİ ARRAYLERİ 1 ARTTIRIYOR
        
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

def select_random_motifs(DNA_lines,k):#RANDOM MOTİF SEÇME / LİST ŞEKLİNDE RETURN EDIYOR
    rows, cols = (len(DNA_lines), k-1)
    motifs  = [[0] * cols] * rows#HER LİNEDAN K-MER KADAR BİR 2D ARRAY OLUŞTURUYOR
    i=0
    
    for lines in DNA_lines:#HER LİNEDAN RANDOM BİR ŞEKİLDE K-MER UZUNLUĞU KADAR STRİNG ALIYOR VE MOTIFS ARRAYINE YAZIYOR
        
        random_position=random.randint(0,len(lines)-k)
        k_mer=lines[random_position:random_position+k]
        motifs[i]=k_mer
        i += 1

    return motifs

def best_match_for_motifs(profile,DNA_lines,k):#
    line=0
    best_mached_motifs_list=([0]*len(DNA_lines))
    best_mached_motifs_score_list = ([0] * len(DNA_lines))
    
    for each_line in DNA_lines:#DNA DAKİ HER LİNE İÇİN FOR A GİRİYOR
        best_match = 0# Her line için farklı best match bulacağımız için sıfırlıyoruz
        
        for y in range(0,len(each_line)-k):#İNPUT LİNE'IN UZUNLUĞU - K-MER KADAR FOR A GİRİYOR
            
            match_for_each_motif = 1#HER MOTİF İÇİN MATCH'E BAKIYORUZ
            motif=each_line[y:y+k]#k-mer şeklinde diziyor

            for motif_index in range(0,k):
                if motif[motif_index]=="A":
                    match_for_each_motif=profile["A"][motif_index]*match_for_each_motif

                    if profile["A"][motif_index]==0:
                        match_for_each_motif=-1
                        break

                if motif[motif_index]=="T":
                    match_for_each_motif=profile["T"][motif_index]*match_for_each_motif
                    if profile["T"][motif_index]==0:
                        match_for_each_motif = -1
                        break

                if motif[motif_index]=="C":
                    match_for_each_motif=profile["C"][motif_index]*match_for_each_motif
                    if profile["C"][motif_index]==0:
                        match_for_each_motif = -1
                        break

                if motif[motif_index]=="G":
                    match_for_each_motif=profile["G"][motif_index]*match_for_each_motif
                    if profile["G"][motif_index]==0:
                        match_for_each_motif = -1
                        break
                    
            if match_for_each_motif>best_match and match_for_each_motif >0:
                best_match=match_for_each_motif
                best_mached_motifs_list[line]=motif
                best_mached_motifs_score_list[line]=best_match
        line+=1

    return best_mached_motifs_list

def randomized_motif_search(DNA,k):

    DNA_lines = DNA.readlines()
    motifs=select_random_motifs(DNA_lines,k)
    bestMotifs = motifs.copy()
    iteration = 0
    while True:
        profil = create_profil(motifs, k,DNA_lines)
        motifs = best_match_for_motifs(profil,DNA_lines,k)
        if find_score(motifs,k,DNA_lines) < find_score(bestMotifs,k,DNA_lines):
            bestMotifs = motifs.copy()
            iteration +=1
        else:
            print("BEST MOTIFS: ",bestMotifs)
            print("SCORE: ",find_score(bestMotifs,k,DNA_lines))
            print("CONSENSUS: ",consensus(bestMotifs,k,DNA_lines))
            print("Iteration number: ",iteration)
            return bestMotifs

def main():
    
    DNA = open("input.txt", "r")
    k = 10

    start = timeit.default_timer()
    randomized_motif_search(DNA,k)
    stop = timeit.default_timer()

    print('Time: ', stop - start)

main()