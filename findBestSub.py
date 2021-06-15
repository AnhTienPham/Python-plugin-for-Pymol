
from pymol import cmd, stored, math
import pandas

## This module identify the best end substitution amino acid for a chosen start amino acid
## based on the abundacne data and charge difference between the two amino acids
## When use in Pymol, must type "import pandas" after download pandas online


## This function returns the score of the amino acid
def aaScore(aa):
    # Acid AA are given -1 score
    # Basic AA are given 1 score
    # Polar AA are given 0 score
    # Hydrophobic AA are given 1 score
    listAcid = ['D','E']
    listBasic = ['R','H','K']
    listPolar = ['S','T','N','Q','C','P']
    listHydro = ['A','V','I','L','M','F','Y','W','G']
    if aa in listAcid:
        return -1
    if aa in listPolar:
        return 0
    if aa in listBasic or listHydro:
        return 1

## Import the abundance data
df = pandas.read_csv('AbundanceData.csv')

## This function identifies the end amino acid with good abundacne score given a start one
def bestTolScore(aa):
    ## Here we choose good abundacne score as > 50, but it can be changed
    cutOff = 50
    listTol = {"End": [], "Score": [], "RankScore": []}
    ## Look into the dataframe which end amino acid has abundance higher than cutoff point
    for ind in df.index:
        if df['start'][ind] == aa and df['pct_tolerated'][ind] > cutOff:
            listTol["End"].append(df['end'][ind])
            listTol["Score"].append(df['pct_tolerated'][ind])
            listTol["RankScore"].append(0)
    listRank = pandas.DataFrame(data=listTol)
    ## Sort the abundance data from lowest to highest
    listRank.sort_values(by=['Score'], inplace=True)
    listRank = listRank.reset_index(drop=True)
    ## Assign a score for each end amino acid based on their abundance ranking
    count = 1
    for ind in listRank.index:
        listRank["RankScore"][ind] = count
        count = count + 1
    return listRank


## This function combines the previous two functions to find the best end amino acids
def findBestSub(aa):
    ## Create a list for the charge difference and get the abundance rank list from the
    ## above function
    listCharge = {"End": [], "chargeDiff": []}
    startAA = bestTolScore(aa)
    ## Calculate the charge difference score for all start-end amino acids pair
    for ind in df.index:
        if df["start"][ind] == aa:
            listCharge["End"].append(df["end"][ind])
            listCharge["chargeDiff"].append(abs(aaScore(aa) - aaScore(df["end"][ind])))
    listCharge = pandas.DataFrame(data=listCharge)
    ## Merge the charge difference and abundance dataframe list into one
    finalList = pandas.merge(listCharge, startAA, on='End')
    ## Create a new column called 'finalScore" that adds up the charge difference
    ## and abundance score to determine the best substitution
    finalList.insert(4, "finalScore", list(range(0, len(finalList))), True)
    for ind in finalList.index:
        finalList["finalScore"][ind] = finalList["chargeDiff"][ind]*50/2 + finalList["RankScore"][ind]*50/len(finalList)
    finalList.sort_values(by=['finalScore'], inplace=True)
    finalList = finalList.reset_index(drop=True)
    ## Here we choose the best 3 substitution, but it can be changed
    numSub = 5
    chosenSubList = []
    ## We choose the best subsitution by taking from the last index of the dataframe list
    ## and print it out
    for x in range(1, numSub + 1):
        chosenSubList.append(finalList["End"][len(finalList) - x])
    print(chosenSubList)

cmd.extend("findBestSub", findBestSub)