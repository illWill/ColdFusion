adapters = [
'AATGATACGGCGACCACCGAGATCTACACGACATTGTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
'AATGATACGGCGACCACCGAGATCTACACACTGATGGACACTCTTTCCCTACACGACGCTCTTCCGATCT',
'AATGATACGGCGACCACCGAGATCTACACGTACCTAGACACTCTTTCCCTACACGACGCTCTTCCGATCT',
'AATGATACGGCGACCACCGAGATCTACACCAGAGCTAACACTCTTTCCCTACACGACGCTCTTCCGATCT',
'AATGATACGGCGACCACCGAGATCTACACCATAGTGAACACTCTTTCCCTACACGACGCTCTTCCGATCT',
'AATGATACGGCGACCACCGAGATCTACACTACCTAGTACACTCTTTCCCTACACGACGCTCTTCCGATCT',
'AATGATACGGCGACCACCGAGATCTACACCGCGATATACACTCTTTCCCTACACGACGCTCTTCCGATCT']
adapter = 'AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT'
water(read.sequence,adapter)

test = [alignment.water(read.sequence[1:150],adapter) for read in reads0.readList]