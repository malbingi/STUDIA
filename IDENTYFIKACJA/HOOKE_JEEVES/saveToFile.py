def saveToFile(function, value):
    dict = {
        0: "Test_func.txt",
        1: "Rosenbrock.txt",
        2: "Rastring.txt"
    }
    
    file = open(dict.get(function), 'a')
    file.write(str(value) + " \n")
    file.close()