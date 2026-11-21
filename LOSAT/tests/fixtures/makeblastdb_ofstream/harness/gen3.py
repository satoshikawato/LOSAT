out=[]
for mode in (0,1):
    for lim in (0,1,5,100,8190,8191,8192,9000,16000):
        for b,n in ((0,10),(0,9000),(100,20),(8000,500),(8190,1),(5000,12000)):
            out += [f"R {lim} {mode}","O",f"W {b}",f"W {n}","F","W 1","F","W 3","F","C","W 1","F","Z"]
open("seq3.txt","w").write("\n".join(out)+"\n")
