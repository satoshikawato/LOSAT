# the state table: every stream state s (0..7) before each operation, on an open stream with
# 5 bytes held (limit -1 and limit 0), on a closed stream, and on a stream with no file
out=[]
opsl=["W 3","W 9000","I 3","P","F","E","S 0","C","O","X","T 1","T 2","T 4"]
for lim,mode in ((-1,0),(0,0)):
    for s in range(8):
        for op in opsl:
            out+=[f"R {lim} {mode}","O","W 5","X",f"T {s}",op,"F","Z"]
for s in range(8):
    for op in opsl:
        out+=["R -1 0","O","W 5","C","X",f"T {s}",op,"Z"]
for s in range(8):
    for op in opsl:
        out+=["R -1 0","OM",f"T {s}",op,"Z"]
open("seq_states.txt","w").write("\n".join(out)+"\n")
