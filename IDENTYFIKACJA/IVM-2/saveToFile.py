def saveToFile(dyslokacje, naprezenia):
    file = open("test1.csv", 'a')
    file.write(str(dyslokacje) + ", " + str(naprezenia) + " \n")
    file.close()