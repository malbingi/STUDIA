def saveToFile(function, value):
    dict = {
        0: "Test_func_rand.txt",
        1: "Rosenbrock_rand.txt",
        2: "Rastring_rand.txt"
    }
    
    file = open(dict.get(function), 'a')
    file.write(str(value) + " \n")
    file.close()