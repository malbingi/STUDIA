def saveToFile(value):
    
    file = open("5.txt", 'a')
    file.write(str(value) + " \n")
    file.close()