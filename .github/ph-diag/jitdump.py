# gdb python: for every thread with a frame inside JSC::JITWorklist::waitUntilAllPlansForVMAreReady
# (libWPEWebKit 2311 vaddr 0x386b650..0x386b91a), print the JITWorklist (r15 of that frame) and its plans.
import gdb, re
LIB = 'libWPEWebKit-2.0.so.1'
base = None
pid = gdb.selected_inferior().pid
for line in open(f'/proc/{pid}/maps'):
    f = line.split()
    if len(f) >= 6 and f[5].endswith(LIB + '.10.2') and int(f[2], 16) == 0:
        base = int(f[0].split('-')[0], 16); break
print(f'jitdump: libWPEWebKit base {base:#x}' if base else 'jitdump: no base')
def u32(a): return int(gdb.parse_and_eval(f'*(unsigned int*){a:#x}'))
def u64(a): return int(gdb.parse_and_eval(f'*(unsigned long*){a:#x}'))
done = set()
for th in gdb.selected_inferior().threads():
    th.switch()
    f = gdb.newest_frame(); n = 0
    while f is not None and n < 40:
        pc = f.pc()
        if base and 0x386b650 <= pc - base < 0x386b91a:
            f.select()
            w = int(f.read_register('r15')) & 0xffffffffffffffff
            vm = int(f.read_register('r14')) & 0xffffffffffffffff
            print(f'jitdump: thread {th.num} lwp {th.ptid[1]} frame {n} worklist {w:#x} vm {vm:#x}')
            if w in done: break
            done.add(w)
            print('  activeThreads', u32(w), 'ongoing', [u32(w + 8 + 4 * i) for i in range(3)], 'max', [u32(w + 0x14 + 4 * i) for i in range(3)])
            print('  raw', ' '.join(f'{u64(w + 8 * i):#x}' for i in range(0x30)))
            for q, off in enumerate((0x40, 0x68, 0x90)):
                print(f'  queue{q} start {u64(w + off):#x} end {u64(w + off + 8):#x} cap/buf {u64(w + off + 0x10):#x} {u64(w + off + 0x18):#x}')
            t = u64(w + 0xb8)
            if t:
                size = u32(t - 4); print('  plans table', hex(t), 'size', size, 'keycount', u32(t - 12), 'deleted', u32(t - 8))
                for i in range(min(size, 4096)):
                    p = u64(t + i * 0x18 + 0x10)
                    if p in (0, 0xffffffffffffffff): continue
                    print(f'   plan {p:#x} key {u64(t + i * 0x18):#x} {u64(t + i * 0x18 + 8):#x} stage {u32(p + 0xc)} vm {u64(p + 0x20):#x} raw', ' '.join(f'{u64(p + 8 * k):#x}' for k in range(8)))
            break
        f = f.older(); n += 1
