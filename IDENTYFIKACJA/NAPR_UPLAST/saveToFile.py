def saveToFile(value, name_id):
    dict  = {
        0: 'model.txt',
        1: '850_01.txt',
        2: '850_1.txt',
        3: '850_10.txt',
        4: '1000_01.txt',
        5: '1000_1.txt',
        6: '1000_10.txt',
        7: '1150_01.txt',
        8: '1150_1.txt',
        9: '1150_10.txt'
    }
    file = open(dict[name_id], 'a')
    file.write(str(value) + " \n")
    file.close()