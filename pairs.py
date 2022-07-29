import string
import itertools

def all_pairs(lst):
    for p in itertools.permutations(lst):
        i = iter(p)
        yield zip(i,i)

def print_states(states, N):
    nstates = 0
    count = 0
    line = ' '
    for s in states:
        l = list(s)
        nstates += 1
        s = '['
        first = True
        for li in l:
            if first == False:
                s += ', '
            s += str(li[0]) + '-' + str(li[1])
            first = False
        s += ']'
        count = count+1

        line += ' ' + s
        if count == N:
            print(line)
            count = 0
            line = ' '
    return nstates