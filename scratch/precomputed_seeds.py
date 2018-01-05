from libnano.seqsearch import seedfinder

START_LIM = 9
END_LIM = 11
print("****** FINDING SINGLES")
for i in range(START_LIM, END_LIM):
    for j in range(1, 3):
        print((i, j), seedfinder.findSeed(i, j))

print("****** FINDING PAIRS")
for i in range(START_LIM, END_LIM):
    for j in range(1, 3):
        print((i, j), seedfinder.findSeedPairs(i, j))

