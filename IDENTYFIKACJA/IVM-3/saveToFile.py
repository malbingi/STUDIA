def saveToFile(dyslokacje, name):
    file = open(name, 'a')
    file.write(str(dyslokacje) + " \n")
    file.close()