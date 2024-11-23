import numpy as np
values_of_K = np.concatenate([np.arange(40,511,8),np.arange(512,1023,16),np.arange(1024,2049,32)])
print(values_of_K)
CRC = 24
Z = 2048


def find_K(bits):
    i = values_of_K.searchsorted(bits)
    if i > 0:
        return values_of_K[i],values_of_K[i-1]
    else:
        return values_of_K[i],None


def segmentation(tbs):
    num_CB = tbs / Z
    
    if num_CB <= 1.0:
#        print("One code block is enough!")
        num_CB = 1        
        k, _ = find_K(tbs)        
        filler = tbs - k
        if filler>0:
            print(f"{filler} bits left!")
        
        return [int(k)]
    
    num_CB = int(np.ceil(num_CB))
#    print("Need > 1 code block")
    
    bits_per_CB = tbs / num_CB
    kplus, kminus = find_K(bits_per_CB)
    
    #print(f"{tbs=} desired {bits_per_CB=} found {kplus=} {kminus=}")
    for num_kplus in range(num_CB, 0,-1):
        left = tbs - num_kplus * kplus
        #print(num_kplus, left, left % kminus)
        if left % kminus == 0:
            break
    else:
        raise RuntimeError("OMG")
    num_kminus = num_CB - num_kplus
    cb = [int(kplus) for i in range(num_kplus)]
    cb += [int(kminus) for i in range(num_kminus)]
    # cb = [int(kplus) for i in range(num_CB-1)]
    # total_bits = sum(cb)
    # left_bits = tbs - total_bits
    # last_k, _ = find_K(left_bits)
    # cb.append(int(last_k))
    
    filler = tbs - sum(cb)
    if filler>0:
        print(f"{filler} bits left!")
    return cb
    
def all_possible_tbs():
    tbs = values_of_K[0]
    while tbs < 40000:
        if tbs < 512:
            tbs += 8
        elif tbs < 1024:
            tbs += 16
        elif tbs < 2048:
            tbs += 32
        else:
            tbs += 64
        yield int(tbs)
            

def desegmentation(tbs):
    num_CB = tbs / Z
    if num_CB < 1.0:
        print("One code block is enough!")
        num_CB = 1        
        k = find_K(tbs)        
        filler = tbs - k
        if filler>0:
            print(f"{filler} bits left!")
        return [int(k)]
    
import csv
csvfile = open('tbs.csv', 'w', newline='')

fieldnames = ['tbs', 'seg1 bits', 'seg1 count', 'seg2 bits', 'seg2 count', 'deviation percent']
writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
writer.writeheader()

#[64*65]:
for tbs in all_possible_tbs():
    segs = segmentation(tbs)
    
    if len(segs)==1:
        usegs = segs
        badness = 0
    else:
        usegs = np.unique(segs)
        
            
    if len(usegs) ==1:
        writer.writerow({"tbs":tbs, "seg1 bits":usegs[0], "seg1 count":len(segs),"seg2 bits":0, "seg2 count":0,"deviation percent":0})
        
    else:
        s1 = max(usegs)
        s2 = min(usegs)
        deviation = int(np.abs(s1 - s2)/s1  * 100)        
        
        n1 = (segs == s1).sum()
        n2 = (segs == s2).sum()
        assert deviation < 5  
        print(f"{tbs=},{n1=},{n2=},{segs=}")
        writer.writerow({"tbs":tbs, "seg1 bits":s1, "seg1 count":n1,"seg2 bits":s2, "seg2 count":n2,"deviation percent":deviation})


csvfile.close()

        
    
