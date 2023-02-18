class Node:

    #VALUES
    score = 0
    match = 3
    gapOpeningPenalty = -1
    gapExtensionPenalty = -0.5
    misMatchPenalty = -1
    isOpeningGap = False

    def __init__(self, xCor, yCor, nucleotid):
        self.xCor = xCor
        self.yCor = yCor
        self.nucleotid = nucleotid
        self.topEdge = 0
        self.diagonalEdge = 0
        self.leftEdge = 0
        
        #Initial nucleotid values.
        if (self.xCor == 0 and self.yCor == 0) or self.xCor!=0 and self.yCor!=0:
            self.nucleotid = "X"

        else:
            self.nucleotid = nucleotid

    #Function for score calculation of each node by choosing maximum score which comes from top, diagonal and left nodes, their edges.
    def calculate_score(self,grid):

        #Score of inital node is 0
        if self.xCor == 0 and self.yCor == 0: 
            self.score = 0

        #Score of first row is increasing -1 by -1.
        elif self.xCor == 0:
            top_x = self.xCor
            top_y = self.yCor - 1
            self.score = grid[top_x][top_y].score + self.gapOpeningPenalty
            self.topEdge = self.gapOpeningPenalty

        #Score of first column is increasing -1 by -1.
        elif self.yCor == 0:
            left_x = self.xCor -1
            left_y = self.yCor
            self.score = grid[left_x][left_y].score + self.gapOpeningPenalty
            self.leftEdge = self.gapOpeningPenalty

        #Calculating scores of nodes except initial node, first row and column.
        else:
            #Getting coordinates of left, top and diagonal nodes.
            left_x = self.xCor - 1; left_y = self.yCor
            top_x = self.xCor; top_y = self.yCor - 1
            diagonal_x = self.xCor - 1; diagonal_y = self.yCor - 1

            max_left_x = 0; max_left_y = self.yCor
            max_x = self.xCor; max_y = 0

            temp_score=[]

            #CHECKING MATCH!
            if(grid[max_left_x][max_left_y].nucleotid==grid[max_x][max_y].nucleotid):
                self.diagonalEdge = self.match
                temp_score.append(grid[diagonal_x][diagonal_y].score+ self.match)

            #CHECKING MISMATCH!
            else:
                temp_score.append(grid[diagonal_x][diagonal_y].score + self.misMatchPenalty)
                self.diagonalEdge = self.misMatchPenalty

            #assigning penalties to the edges which connects nodes.
            self.topEdge = self.gapOpeningPenalty
            self.leftEdge = self.gapOpeningPenalty

            #Selecting maximum score to put the score value according to dynamic programming rules.
            temp_score.append(grid[left_x][left_y].score+self.leftEdge)#Soldan indel
            temp_score.append(grid[top_x][top_y].score+self.topEdge)#Üstten indel
            self.score=max(temp_score)

#backtracking function which backtracks from end of the grid to starting point of grid according to dynamic pairwise alignment rules.
def backTracking(grid,seq1,seq2):

    x = len(seq1)
    y = len(seq2)

    seq1_align = ""
    seq2_align = ""

    gap_extension_scores = 0 # this extension scores will be added to total scores.

    while x+y != 0:
        
        '''
        Main idea of backtracking with gap openings and extension scores is that preferring indel penalties to mismatches in
        specific situation such as putting indel and making mismatch are available at same point because
        gap extension penalty is less than mismatch, so if there will be more than one indel, penalty will be less than mismatches.
        '''
        
        #Checking top edge of the node (INDEL)
        if y != 0 and grid[x][y].score - grid[x][y].topEdge == grid[x][y-1].score:
            seq1_align += '-'
            seq2_align += grid[0][y].nucleotid

            #gap opened
            grid[x][y].isOpeningGap = True
            
            #checking if gap opened before, if yes than it is gap extension.
            if y != len(seq2):
                if grid[x][y+1].isOpeningGap == True:
                    gap_extension_scores += 0.5
            y -= 1

        #Checking left edge of the node (INDEL)
        elif x != 0 and grid[x][y].score - grid[x][y].leftEdge == grid[x-1][y].score: 
            seq1_align += grid[x][0].nucleotid
            seq2_align += '-'

            #gap opened
            grid[x][y].isOpeningGap = True

            #checking if gap opened before, if yes than it is gap extension.
            if x != len(seq1):
                if grid[x+1][y].isOpeningGap == True:
                    gap_extension_scores += 0.5
            x -= 1
        
        #Checking if there is a mismatch or match.
        elif (x != 0 or y != 0) and grid[x][y].score - grid[x][y].diagonalEdge == grid[x-1][y-1].score: #çapraz
            seq1_align += grid[x][0].nucleotid
            seq2_align += grid[0][y].nucleotid
            x -= 1
            y -= 1

    #print the results
    print(seq1_align[::-1])
    print(seq2_align[::-1])
    print("Score:", grid[-1][-1].score + gap_extension_scores)
        
def main():
    test_inputs=["test1.seq","test2.seq","test3.seq","test4.seq","test5.seq"]

    for file in test_inputs:
        with open("Test Inputs/"+file, "r") as f:
            seq1 = f.readline()[:-1]
            seq2 = f.readline()[:-1]
        
        grid = [ [Node] * (len(seq2)+1) for i in range(len(seq1)+1)]
        grid[0][0] = Node(0, 0, 'X')

        #Calculating initial score of second sequence (-1 indels)
        for i in range(1,len(seq2)+1):
            node = Node(0,i,seq2[i-1])
            grid[0][i] = node
            node.calculate_score(grid)

        #Calculating initial score of first sequence (-1 indels)
        for i in range(1,len(seq1)+1):
            node = Node(i,0,seq1[i-1])
            grid[i][0] = node
            node.calculate_score(grid)

        #putting other nodes to the grid and calculates the scores
        for i in range(1,len(seq1)+1):
            for y in range(1, len(seq2) + 1):
                node = Node(i, y, "X")
                grid[i][y] = node
                node.calculate_score(grid)

        #printing the results
        print("")
        print("Result of "+ file)
        backTracking(grid,seq1,seq2)
        
        
main()
